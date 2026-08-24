#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

#include <zspace/zToolsets/fabrication/z3DPSlicer.h>

namespace zSpace
{
    namespace
    {
        class z3dpTimer
        {
        public:
            z3dpTimer()
                : m_start(std::chrono::high_resolution_clock::now())
            {
            }

            double elapsedMs() const
            {
                const auto now = std::chrono::high_resolution_clock::now();
                return std::chrono::duration<double, std::milli>(now - m_start).count();
            }

        private:
            std::chrono::high_resolution_clock::time_point m_start;
        };

        struct z3dpLoopSampleData
        {
            zPointArray positions;
            zFloatArray cumulativeLengths;
            float totalLength = 0.0f;
        };

        using z3dpHalfedgeVectorCache = std::unordered_map<int, zVector>;

        float z3dpMax3(float a, float b, float c)
        {
            return std::max(a, std::max(b, c));
        }

        bool z3dpTryNext(zItMeshHalfEdge &he, zItMeshHalfEdge& out)
        {
            try
            {
                out = he.getNext();
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTryNextSymNext(zItMeshHalfEdge &he, zItMeshHalfEdge& out)
        {
            try
            {
                zItMeshHalfEdge heNext = he.getNext();
                zItMeshHalfEdge heSym = heNext.getSym();
                out = heSym.getNext();
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTrySymPrevSym(zItMeshHalfEdge &he, zItMeshHalfEdge& out)
        {
            try
            {
                zItMeshHalfEdge heSym = he.getSym();
                zItMeshHalfEdge hePrev = heSym.getPrev();
                out = hePrev.getSym();
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTrySym(zItMeshHalfEdge &he, zItMeshHalfEdge& out)
        {
            try
            {
                out = he.getSym();
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTrySymNext(zItMeshHalfEdge &he, zItMeshHalfEdge& out)
        {
            try
            {
                zItMeshHalfEdge heSym = he.getSym();
                out = heSym.getNext();
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTryColorHalfedge(zItMeshHalfEdge &he, const zColor& color, double weight)
        {
            try
            {
                he.getEdge().setColor(color);
                he.getEdge().setWeight(weight);
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTryEndVertexId(zItMeshHalfEdge &he, int& vertexId)
        {
            try
            {
                vertexId = he.getVertex().getId();
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTryEndVertexValence(zItMeshHalfEdge &he, int expectedValence)
        {
            try
            {
                return he.getVertex().getValence() == expectedValence;
            }
            catch (...)
            {
                return false;
            }
        }

        bool z3dpTryAdvanceBlueLoop(zItMeshHalfEdge& he, zItMeshHalfEdge& out)
        {
            const bool isCorner = z3dpTryEndVertexValence(he, 2);
            return isCorner ? z3dpTryNext(he, out) : z3dpTryNextSymNext(he, out);
        }

        bool z3dpTryCachedHalfedgeVector(zItMeshHalfEdge& he, z3dpHalfedgeVectorCache& vectorCache, zVector& out)
        {
            try
            {
                const int heId = he.getId();
                const auto it = vectorCache.find(heId);
                if (it != vectorCache.end())
                {
                    out = it->second;
                    return true;
                }

                out = he.getVector();
                vectorCache[heId] = out;
                return true;
            }
            catch (...)
            {
                return false;
            }
        }

        void z3dpCreateRectangleGraph(zObjectGraph& graph, const zPoint& minBounds, const zPoint& maxBounds, float offset)
        {
            const float maxOffset = std::min(std::abs(maxBounds.x - minBounds.x), std::abs(maxBounds.y - minBounds.y)) * 0.45f;
            const float safeOffset = std::min(std::max(offset, 0.0f), maxOffset);
            zPointArray positions = {
                zPoint(minBounds.x + safeOffset, minBounds.y + safeOffset, 0.0f),
                zPoint(maxBounds.x - safeOffset, minBounds.y + safeOffset, 0.0f),
                zPoint(maxBounds.x - safeOffset, maxBounds.y - safeOffset, 0.0f),
                zPoint(minBounds.x + safeOffset, maxBounds.y - safeOffset, 0.0f)
            };
            zIntArray connects = { 0, 1, 1, 2, 2, 3, 3, 0 };
            zFnGraph fn(graph);
            fn.clear();
            fn.create(positions, connects);
        }

        zObjectGraph z3dpGraphFromHalfedgeLoop(const zItMeshHalfEdgeArray& loop)
        {
            zObjectGraph graph;
            zPointArray positions;
            zIntArray connects;

            for (int i = 0; i < static_cast<int>(loop.size()); i++)
            {
                zItMeshHalfEdge loopHe = loop[i];
                try
                {
                    if (i == 0) positions.push_back(loopHe.getVertex().getPosition());
                    positions.push_back(loopHe.getStartVertex().getPosition());
                }
                catch (...)
                {
                    return graph;
                }
            }

            if (positions.size() < 2) return graph;

            connects.reserve((positions.size() - 1) * 2);
            for (int i = 1; i < static_cast<int>(positions.size()); i++)
            {
                connects.push_back(i - 1);
                connects.push_back(i);
            }

            zFnGraph fn(graph);
            fn.create(positions, connects);
            return graph;
        }

        bool z3dpFaceVerticesFromHalfedge(zItMeshHalfEdge heStart, bool forward, zIntArray& faceVertexIds);

        bool z3dpBuildMappedFacesFromSourcePolygons(
            zObjectMesh& sourceMesh,
            const std::unordered_map<int, int>& sourceVertexToBottomVertex,
            zIntArray& polyCounts,
            zIntArray& polyConnects,
            zIntArray& bottomStripFaceIds)
        {
            zIntArray sourceConnects;
            zIntArray sourceCounts;
            try
            {
                zFnMesh fnSource(sourceMesh);
                fnSource.getPolygonData(sourceConnects, sourceCounts);
            }
            catch (...)
            {
                return false;
            }

            int cursor = 0;
            for (int faceId = 0; faceId < static_cast<int>(sourceCounts.size()); faceId++)
            {
                const int faceCount = sourceCounts[faceId];
                if (faceCount < 3 || cursor + faceCount > static_cast<int>(sourceConnects.size()))
                {
                    cursor += std::max(faceCount, 0);
                    continue;
                }

                zIntArray faceConnects;
                bool validFace = true;
                for (int j = 0; j < faceCount; j++)
                {
                    const int sourceVertexId = sourceConnects[cursor + j];
                    const auto mapIt = sourceVertexToBottomVertex.find(sourceVertexId);
                    if (mapIt == sourceVertexToBottomVertex.end())
                    {
                        validFace = false;
                        break;
                    }
                    faceConnects.push_back(mapIt->second);
                }

                if (validFace)
                {
                    polyCounts.push_back(faceCount);
                    for (const int vertexId : faceConnects) polyConnects.push_back(vertexId);
                    bottomStripFaceIds.push_back(faceId);
                }

                cursor += faceCount;
            }

            return !polyCounts.empty();
        }

        float z3dpPointDistance(const zPoint& a, const zPoint& b)
        {
            const float dx = a.x - b.x;
            const float dy = a.y - b.y;
            const float dz = a.z - b.z;
            return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
        }

        zPoint z3dpLerpPoint(const zPoint& a, const zPoint& b, float t)
        {
            return zPoint(
                a.x + ((b.x - a.x) * t),
                a.y + ((b.y - a.y) * t),
                a.z + ((b.z - a.z) * t));
        }

        float z3dpPolylineLength(const zPointArray& positions)
        {
            float length = 0.0f;
            for (int i = 1; i < static_cast<int>(positions.size()); i++)
            {
                length += z3dpPointDistance(positions[i], positions[i - 1]);
            }
            return length;
        }

        bool z3dpBuildLoopSampleData(zFnMesh& fnMesh, zItMeshHalfEdgeArray& loop, int divisions, bool useStartVertexAsBottom, z3dpLoopSampleData& sampleData)
        {
            sampleData.positions.clear();
            sampleData.cumulativeLengths.clear();
            sampleData.totalLength = 0.0f;

            if (loop.empty() || divisions < 1) return false;

            try
            {
                fnMesh.computeEdgeLoop_Split(loop, divisions, sampleData.positions);
            }
            catch (...)
            {
                return false;
            }

            if (!useStartVertexAsBottom) std::reverse(sampleData.positions.begin(), sampleData.positions.end());

            const zPointArray& positions = sampleData.positions;
            if (positions.size() < 2) return false;

            sampleData.cumulativeLengths.reserve(positions.size());
            sampleData.cumulativeLengths.push_back(0.0f);

            for (int i = 1; i < static_cast<int>(positions.size()); i++)
            {
                sampleData.totalLength += z3dpPointDistance(positions[i], positions[i - 1]);
                sampleData.cumulativeLengths.push_back(sampleData.totalLength);
            }

            return sampleData.totalLength > 1.0e-6f;
        }

        zPoint z3dpSamplePolyline(const z3dpLoopSampleData& sampleData, float t)
        {
            const zPointArray& positions = sampleData.positions;
            if (positions.empty()) return zPoint();
            if (positions.size() == 1) return positions[0];

            t = std::min(std::max(t, 0.0f), 1.0f);

            if (sampleData.totalLength <= 1.0e-6f) return positions[0];

            const float targetLength = sampleData.totalLength * t;

            for (int i = 1; i < static_cast<int>(positions.size()); i++)
            {
                const float previousLength = sampleData.cumulativeLengths[i - 1];
                const float currentLength = sampleData.cumulativeLengths[i];
                const float segmentLength = currentLength - previousLength;
                if (segmentLength <= 1.0e-6f) continue;

                if (currentLength >= targetLength)
                {
                    const float segmentT = (targetLength - previousLength) / segmentLength;
                    return z3dpLerpPoint(positions[i - 1], positions[i], segmentT);
                }
            }

            return positions.back();
        }

        bool z3dpFaceVerticesFromHalfedge(zItMeshHalfEdge heStart, bool forward, zIntArray& faceVertexIds)
        {
            faceVertexIds.clear();
            zItMeshHalfEdge he = heStart;

            try
            {
                int guard = 0;
                do
                {
                    faceVertexIds.push_back(forward ? he.getVertex().getId() : he.getStartVertex().getId());
                    he = forward ? he.getNext() : he.getPrev();
                } while (he.getId() != heStart.getId() && guard++ < 64);
            }
            catch (...)
            {
                return false;
            }

            return faceVertexIds.size() >= 3;
        }

        bool z3dpCreateBottomMeshFromBlendLogic(
            zObjectMesh& sourceMesh,
            const std::vector<zItMeshHalfEdgeArray>& vLoops,
            zItMeshHalfEdgeArray& bottomLoop,
            zObjectMesh& bottomMesh,
            zIntArray& bottomStripFaceIds,
            bool& useStartVertexAsBottom)
        {
            bottomMesh = zObjectMesh();
            bottomStripFaceIds.clear();
            useStartVertexAsBottom = false;
            if (vLoops.empty() || bottomLoop.empty()) return false;

            for (int endpointPass = 0; endpointPass < 2; endpointPass++)
            {
                const bool useStartEndpoint = endpointPass == 1;

                zPointArray positions;
                zIntArray polyCounts;
                zIntArray polyConnects;
                std::unordered_map<int, int> sourceVertexToBottomVertex;

                positions.reserve(vLoops.size());
                for (int i = 0; i < static_cast<int>(vLoops.size()); i++)
                {
                    if (vLoops[i].empty()) continue;

                    try
                    {
                        zItMeshHalfEdge loopStart = vLoops[i][0];
                        const int sourceVertexId = useStartEndpoint ? loopStart.getStartVertex().getId() : loopStart.getVertex().getId();
                        if (sourceVertexToBottomVertex.find(sourceVertexId) != sourceVertexToBottomVertex.end()) continue;

                        sourceVertexToBottomVertex[sourceVertexId] = static_cast<int>(positions.size());
                        positions.push_back(useStartEndpoint ? loopStart.getStartVertex().getPosition() : loopStart.getVertex().getPosition());
                    }
                    catch (...)
                    {
                        continue;
                    }
                }

                if (positions.empty()) continue;
                if (!z3dpBuildMappedFacesFromSourcePolygons(sourceMesh, sourceVertexToBottomVertex, polyCounts, polyConnects, bottomStripFaceIds)) continue;

                zFnMesh fn(bottomMesh);
                fn.create(positions, polyCounts, polyConnects);
                fn.setFaceColor(zColor(0.1f, 0.45f, 1.0f, 0.18f));
                fn.setEdgeColor(zColor(0.1f, 0.55f, 1.0f, 1.0f));
                fn.setEdgeWeight(1.5);
                useStartVertexAsBottom = useStartEndpoint;
                return true;
            }

            bottomStripFaceIds.clear();
            return false;
        }

        void z3dpBuildSliceMeshes(
            zObjectMesh& sourceMesh,
            std::vector<zItMeshHalfEdgeArray>& vLoops,
            zItMeshHalfEdgeArray& bottomLoop,
            const z3DPSettings& settings,
            int& sliceCount,
            zObjectMeshArray& sliceMeshes,
            z3DPSlicingFeatures& features)
        {
            sliceCount = 0;
            sliceMeshes.clear();
            features.topMesh = zObjectMesh();
            features.bottomMesh = zObjectMesh();
            features.bottomStripFaceIds.clear();
            features.topStripFaceIds.clear();

            if (vLoops.size() < 2) return;

            z3dpTimer totalTimer;

            z3dpTimer bottomTimer;
            zIntArray bottomStripFaceIds;
            zObjectMesh bottomMesh;
            bool useStartVertexAsBottom = false;
            if (!z3dpCreateBottomMeshFromBlendLogic(sourceMesh, vLoops, bottomLoop, bottomMesh, bottomStripFaceIds, useStartVertexAsBottom)) return;
            const double bottomMs = bottomTimer.elapsedMs();

            z3dpTimer loopTimer;
            zFnMesh fnSource(sourceMesh);

            float shortestLoopLength = std::numeric_limits<float>::max();
            for (auto& loop : vLoops)
            {
                if (loop.empty()) continue;

                try
                {
                    const float loopLength = fnSource.computeEdgeLoop_Length(loop);
                    if (loopLength > 1.0e-6f) shortestLoopLength = std::min(shortestLoopLength, loopLength);
                }
                catch (...)
                {
                    continue;
                }
            }
            const double loopMs = loopTimer.elapsedMs();

            if (shortestLoopLength == std::numeric_limits<float>::max()) return;

            const float layerHeight = std::max(settings.printLayerHeight, 0.0001f);
            sliceCount = std::max(1, static_cast<int>(std::ceil(shortestLoopLength / layerHeight)));
            const int loopDivisions = std::max(1, sliceCount - 1);

            z3dpTimer splitTimer;
            std::vector<zPointArray> dividedLoopPoints;
            dividedLoopPoints.reserve(vLoops.size());
            for (auto& loop : vLoops)
            {
                z3dpLoopSampleData sampleData;
                if (!z3dpBuildLoopSampleData(fnSource, loop, loopDivisions, useStartVertexAsBottom, sampleData)) continue;
                dividedLoopPoints.push_back(sampleData.positions);
            }
            const double splitMs = splitTimer.elapsedMs();

            if (dividedLoopPoints.size() < 2) return;
            sliceCount = static_cast<int>(dividedLoopPoints.front().size());
            sliceMeshes.reserve(sliceCount);

            z3dpTimer topologyTimer;
            zFnMesh fnBottom(bottomMesh);
            zIntArray bottomConnects;
            zIntArray bottomCounts;
            fnBottom.getPolygonData(bottomConnects, bottomCounts);
            const double topologyMs = topologyTimer.elapsedMs();

            double sampleMs = 0.0;
            double createMs = 0.0;
            double styleMs = 0.0;
            double pushMs = 0.0;

            for (int i = 0; i < sliceCount; i++)
            {
                z3dpTimer sampleTimer;
                zPointArray positions;
                positions.reserve(dividedLoopPoints.size());
                for (const auto& loopPoints : dividedLoopPoints)
                {
                    if (i < static_cast<int>(loopPoints.size())) positions.push_back(loopPoints[i]);
                }
                sampleMs += sampleTimer.elapsedMs();

                z3dpTimer createTimer;
                zObjectMesh sliceMesh;
                zFnMesh fnSlice(sliceMesh);
                fnSlice.create(positions, bottomCounts, bottomConnects);
                createMs += createTimer.elapsedMs();

                z3dpTimer styleTimer;
                fnSlice.setFaceColor(zColor(0.1f, 0.45f, 1.0f, 0.18f));
                fnSlice.setEdgeColor(zColor(0.1f, 0.55f, 1.0f, 1.0f));
                fnSlice.setEdgeWeight(1.5);
                styleMs += styleTimer.elapsedMs();

                z3dpTimer pushTimer;
                sliceMeshes.push_back(sliceMesh);
                pushMs += pushTimer.elapsedMs();
            }

            if (!sliceMeshes.empty())
            {
                features.bottomMesh = sliceMeshes.front();
                features.topMesh = sliceMeshes.back();
                features.bottomStripFaceIds = bottomStripFaceIds;
                features.topStripFaceIds = bottomStripFaceIds;
            }

            std::cout << "[3DP timing] computeSlices detail"
                << " | bottom topology: " << bottomMs << " ms"
                << " | loop lengths: " << loopMs << " ms"
                << " | split loop points: " << splitMs << " ms"
                << " | polygon data: " << topologyMs << " ms"
                << " | sample vertices: " << sampleMs << " ms"
                << " | create meshes: " << createMs << " ms"
                << " | style meshes: " << styleMs << " ms"
                << " | store meshes: " << pushMs << " ms"
                << " | slices: " << sliceCount
                << " | loops: " << dividedLoopPoints.size()
                << " | total: " << totalTimer.elapsedMs() << " ms"
                << std::endl;
        }

        bool z3dpWalkVerticalLoop(zItMeshHalfEdge heStart, int verticalStride, zItMeshHalfEdgeArray& loop)
        {
            loop.clear();
            if (verticalStride <= 0) return false;

            zItMeshHalfEdge he = heStart;
            for (int i = 0; i < verticalStride; i++)
            {
                he.getEdge().setColor(zColor(1, 0, 1, 1));
                he.getEdge().setWeight(4.0);
                loop.push_back(he);

                if (i == verticalStride - 1) break;

                he = he.getNext().getSym().getNext();
            }

            return !loop.empty();
        }

        bool z3dpTryBestVerticalStart(zItMeshHalfEdge& blueHe, const zVector& verticalDir, z3dpHalfedgeVectorCache& vectorCache, zItMeshHalfEdge& out)
        {
            zItMeshHalfEdgeArray candidates;
            zVector referenceDir = verticalDir;

            zItMeshHalfEdge candidate;
            if (z3dpTrySymPrevSym(blueHe, candidate)) candidates.push_back(candidate);
            if (z3dpTryNext(blueHe, candidate)) candidates.push_back(candidate);
            if (z3dpTrySymNext(blueHe, candidate)) candidates.push_back(candidate);

            const int candidateCount = static_cast<int>(candidates.size());
            for (int i = 0; i < candidateCount; i++)
            {
                zItMeshHalfEdge candidateSym;
                if (z3dpTrySym(candidates[i], candidateSym)) candidates.push_back(candidateSym);
            }

            bool found = false;
            float bestAngle = std::numeric_limits<float>::max();

            for (auto& he : candidates)
            {
                zVector heDir;
                if (!z3dpTryCachedHalfedgeVector(he, vectorCache, heDir)) continue;
                if (heDir.length() <= 1.0e-6f) continue;

                const float angle = heDir.angle(referenceDir);
                if (angle < bestAngle)
                {
                    bestAngle = angle;
                    out = he;
                    found = true;
                }
            }

            return found;
        }

        bool z3dpWalkBlueLoop(zItMeshHalfEdge heStart, int startVID, int horizontalVID, int maxSteps, zItMeshHalfEdgeArray& loop, bool colorEdges)
        {
            loop.clear();
            bool containsHorizontalSeed = false;
            zItMeshHalfEdge he = heStart;

            for (int i = 0; i < maxSteps; i++)
            {
                if (colorEdges && !z3dpTryColorHalfedge(he, zColor(0, 0, 1, 1), 4.0)) break;
                loop.push_back(he);

                int endVID = -1;
                if (!z3dpTryEndVertexId(he, endVID)) break;
                if (endVID == horizontalVID) containsHorizontalSeed = true;
                if (i > 0 && endVID == startVID) break;

                zItMeshHalfEdge heNext;
                if (!z3dpTryAdvanceBlueLoop(he, heNext)) break;
                he = heNext;
            }

            return !loop.empty() && containsHorizontalSeed;
        }

        void z3dpColorHalfedgeLoop(zItMeshHalfEdgeArray& loop, const zColor& color, double weight)
        {
            for (auto& he : loop)
            {
                z3dpTryColorHalfedge(he, color, weight);
            }
        }

        bool z3dpTryBestBlueLoopStart(zItMeshHalfEdgeArray& candidates, const zVector& horizontalDir, int startVID, int horizontalVID, int maxSteps, z3dpHalfedgeVectorCache& vectorCache, zItMeshHalfEdge& outStart, zItMeshHalfEdgeArray& outLoop)
        {
            bool found = false;
            bool foundSeededLoop = false;
            float bestScore = std::numeric_limits<float>::max();
            zVector referenceDir = horizontalDir;

            for (auto& candidate : candidates)
            {
                zVector candidateDir;
                if (!z3dpTryCachedHalfedgeVector(candidate, vectorCache, candidateDir)) continue;
                if (candidateDir.length() <= 1.0e-6f) continue;

                zItMeshHalfEdgeArray candidateLoop;
                const bool hasHorizontalSeed = z3dpWalkBlueLoop(candidate, startVID, horizontalVID, maxSteps, candidateLoop, false);

                const float angle = candidateDir.angle(referenceDir);
                const float score = angle + (hasHorizontalSeed ? 0.0f : 10000.0f);
                if (score < bestScore)
                {
                    bestScore = score;
                    outStart = candidate;
                    outLoop = candidateLoop;
                    found = true;
                    foundSeededLoop = hasHorizontalSeed;
                }
            }

            return found && foundSeededLoop && !outLoop.empty();
        }

        void z3dpComputeVLoops(zObjectMesh& oMesh, int startVID, int endVID, int horizontalVID, vector<zItMeshHalfEdgeArray>& v_Loops, zItMeshHalfEdgeArray& bottomLoop)
        {
            z3dpTimer totalTimer;
            double validateMs = 0.0;
            double seedHalfedgesMs = 0.0;
            double firstLoopMs = 0.0;
            double blueLoopMs = 0.0;
            double blueColorMs = 0.0;
            double verticalLoopMs = 0.0;
            double verticalStartMs = 0.0;
            double verticalWalkMs = 0.0;

            z3dpTimer validateTimer;
            bottomLoop.clear();
            zFnMesh fnMesh_in(oMesh);
            const int numV = fnMesh_in.numVertices();
            if (startVID < 0 || startVID >= numV || endVID < 0 || endVID >= numV || horizontalVID < 0 || horizontalVID >= numV) return;
            z3dpHalfedgeVectorCache halfedgeVectorCache;
            halfedgeVectorCache.reserve(std::max(1, fnMesh_in.numEdges() * 2));

            zItMeshVertex vStart(oMesh, startVID);
            zItMeshVertex vEnd(oMesh, endVID);
            zItMeshVertex vHorizontal(oMesh, horizontalVID);

            zVector dir = vEnd.getPosition() - vStart.getPosition();
            zVector dir_horizontal = vHorizontal.getPosition() - vStart.getPosition();
            validateMs = validateTimer.elapsedMs();

            z3dpTimer seedHalfedgesTimer;
            zItMeshHalfEdgeArray hEdges_Start;
            vStart.getConnectedHalfEdges(hEdges_Start);

            float ang = 10000.0f;
            zItMeshHalfEdge heStart;
            bool foundHeStart = false;
            
            for (auto& he : hEdges_Start)
            {
                zVector vec;
                if (!z3dpTryCachedHalfedgeVector(he, halfedgeVectorCache, vec)) continue;
                if (vec.length() > 1.0e-6f)
                {
                    float aVal = vec.angle(dir);
                    if (aVal < ang)
                    {
                        ang = aVal;
                        heStart = he;
                        foundHeStart = true;
                    }
                    
                }
            }

            if (!foundHeStart) return;
            seedHalfedgesMs = seedHalfedgesTimer.elapsedMs();

            z3dpTimer firstLoopTimer;
            zItMeshHalfEdge he_v = heStart;
            zItMeshHalfEdgeArray firstVerticalLoop;
            int steps_v = 0;
            
            while (steps_v++ < 200)
            {
                he_v.getEdge().setColor(zColor(1, 0, 1, 1));
                he_v.getEdge().setWeight(4.0);
                firstVerticalLoop.push_back(he_v);

                int v_end = -1;
                if (!z3dpTryEndVertexId(he_v, v_end)) break;
                if (v_end == endVID) break;

                he_v = he_v.getNext().getSym().getNext();
            }

            const int verticalStride = static_cast<int>(firstVerticalLoop.size());
            if (verticalStride == 0) return;
            v_Loops.push_back(firstVerticalLoop);
            firstLoopMs = firstLoopTimer.elapsedMs();

            z3dpTimer blueLoopTimer;
            hEdges_Start.clear();
            vStart.getConnectedHalfEdges(hEdges_Start);

            zItMeshHalfEdgeArray horizontalLoop;
            zItMeshHalfEdge heStart_horizontal;
            const int maxHorizontalSteps = std::max(1, fnMesh_in.numEdges() + 4);
            if (!z3dpTryBestBlueLoopStart(hEdges_Start, dir_horizontal, startVID, horizontalVID, maxHorizontalSteps, halfedgeVectorCache, heStart_horizontal, horizontalLoop)) return;
            blueLoopMs = blueLoopTimer.elapsedMs();

            z3dpTimer blueColorTimer;
            z3dpColorHalfedgeLoop(horizontalLoop, zColor(0, 0, 1, 1), 4.0);
            bottomLoop = horizontalLoop;
            blueColorMs = blueColorTimer.elapsedMs();

            z3dpTimer verticalLoopTimer;
            std::unordered_set<int> verticalSeedIds;
            verticalSeedIds.insert(firstVerticalLoop.front().getId());

            for (auto& blueHe : horizontalLoop)
            {
                zItMeshHalfEdge verticalStart;
                z3dpTimer verticalStartTimer;
                if (!z3dpTryBestVerticalStart(blueHe, dir, halfedgeVectorCache, verticalStart)) continue;
                verticalStartMs += verticalStartTimer.elapsedMs();

                const int seedId = verticalStart.getId();
                if (verticalSeedIds.find(seedId) != verticalSeedIds.end()) continue;

                zItMeshHalfEdgeArray verticalLoop;
                z3dpTimer verticalWalkTimer;
                if (!z3dpWalkVerticalLoop(verticalStart, verticalStride, verticalLoop)) continue;
                verticalWalkMs += verticalWalkTimer.elapsedMs();

                verticalSeedIds.insert(seedId);
                v_Loops.push_back(verticalLoop);
            }
            verticalLoopMs = verticalLoopTimer.elapsedMs();

            std::cout << "[3DP timing] compute loops detail"
                << " | validate/seeds: " << validateMs << " ms"
                << " | start halfedges: " << seedHalfedgesMs << " ms"
                << " | first vertical: " << firstLoopMs << " ms"
                << " | blue loop search: " << blueLoopMs << " ms"
                << " | blue color: " << blueColorMs << " ms"
                << " | vertical loops: " << verticalLoopMs << " ms"
                << " | vertical starts: " << verticalStartMs << " ms"
                << " | vertical walks: " << verticalWalkMs << " ms"
                << " | loops: " << v_Loops.size()
                << " | bottom loop edges: " << bottomLoop.size()
                << " | total: " << totalTimer.elapsedMs() << " ms"
                << std::endl;

        }
    }

    void z3DPSlicer::setInputMesh(zObjectMesh& mesh)
    {
        m_inputMesh = &mesh;
    }

    void z3DPSlicer::setSettings(const z3DPSettings& settings)
    {
        m_settings = settings;
    }

    void z3DPSlicer::identifySlicingFeatures()
    {
        if (!m_inputMesh) return;
        z3dpTimer totalTimer;

        z3dpTimer boundsTimer;
        zFnMesh fn(*m_inputMesh);
        fn.getBounds(m_minBounds, m_maxBounds);
        const double boundsMs = boundsTimer.elapsedMs();

        z3dpTimer clearTimer;
        m_features = z3DPSlicingFeatures();
        m_sliceMeshes.clear();
        m_triangulatedMeshes.clear();
        m_unrolledMeshes.clear();
        m_fields.clear();
        m_contours.clear();
        m_verticalLoops.clear();
        m_bottomLoop.clear();
        const double clearMs = clearTimer.elapsedMs();

        z3dpTimer edgeResetTimer;
        fn.setEdgeColor(zColor(0.62f, 0.62f, 0.66f, 1.0f));
        fn.setEdgeWeight(1.0);
        const double edgeResetMs = edgeResetTimer.elapsedMs();

        int fromVertexId = m_settings.slicingSeedFromVertexId;
        int toVertexId = m_settings.slicingSeedToVertexId;
        int horizontalVertexId = m_settings.slicingSeedHorizontalVertexId;

        const int numV = fn.numVertices();
        if (fromVertexId < 0 || fromVertexId >= numV || toVertexId < 0 || toVertexId >= numV || horizontalVertexId < 0 || horizontalVertexId >= numV) return;

        // 1. Compute longitudinal columns (vLoops).
        z3dpTimer loopsTimer;
        z3dpComputeVLoops(*m_inputMesh, fromVertexId, toVertexId, horizontalVertexId, m_verticalLoops, m_bottomLoop);
        const double loopsMs = loopsTimer.elapsedMs();

        // 2. Store loop count only. Vertical loops remain as halfedge arrays;
        // visualization reads edge color/weight directly from the input mesh.
        z3dpTimer graphTimer;
        m_features.edgeLoops.clear();
        m_features.edgeLoopCount = static_cast<int>(m_verticalLoops.size());
        const double graphMs = graphTimer.elapsedMs();

        std::cout << "[3DP timing] identifySlicingFeatures detail"
            << " | bounds: " << boundsMs << " ms"
            << " | clear: " << clearMs << " ms"
            << " | edge reset: " << edgeResetMs << " ms"
            << " | compute loops: " << loopsMs << " ms"
            << " | graph objects skipped: " << graphMs << " ms"
            << " | loops: " << m_verticalLoops.size()
            << " | total: " << totalTimer.elapsedMs() << " ms"
            << std::endl;
    }

    void z3DPSlicer::computeSliceCount()
    {
        const float dx = std::abs(m_maxBounds.x - m_minBounds.x);
        const float dy = std::abs(m_maxBounds.y - m_minBounds.y);
        const float dz = std::abs(m_maxBounds.z - m_minBounds.z);
        const float printLength = z3dpMax3(dx, dy, dz);
        const float h = std::max(m_settings.printLayerHeight, 0.0001f);
        m_sliceCount = std::max(1, static_cast<int>(std::ceil(printLength / h)));
    }

    void z3DPSlicer::computeSliceMeshes()
    {
        if (!m_inputMesh) return;
        if (m_verticalLoops.empty() || m_bottomLoop.empty()) identifySlicingFeatures();
        z3dpBuildSliceMeshes(*m_inputMesh, m_verticalLoops, m_bottomLoop, m_settings, m_sliceCount, m_sliceMeshes, m_features);
    }

    void z3DPSlicer::triangulateSlices()
    {
        m_triangulatedMeshes = m_sliceMeshes;
    }

    void z3DPSlicer::unrollSlices()
    {
        m_unrolledMeshes = m_triangulatedMeshes.empty() ? m_sliceMeshes : m_triangulatedMeshes;
        for (auto& mesh : m_unrolledMeshes)
        {
            zFnMesh fn(mesh);
            zPointArray positions;
            fn.getVertexPositions(positions);
            for (auto& p : positions) p.z = 0.0f;
            fn.setVertexPositions(positions);
        }
    }

    void z3DPSlicer::computeSDFFields()
    {
        m_fields.clear();
        m_fields.assign(m_unrolledMeshes.size(), zObjectMeshScalarField());

        for (int i = 0; i < static_cast<int>(m_unrolledMeshes.size()); ++i)
        {
            zFnMesh fnMesh(m_unrolledMeshes[i]);
            zPoint minBounds;
            zPoint maxBounds;
            fnMesh.getBounds(minBounds, maxBounds);

            zFnMeshScalarField fnField(m_fields[i]);
            fnField.create(minBounds, maxBounds, m_settings.fieldResolutionX, m_settings.fieldResolutionY);
        }
    }

    void z3DPSlicer::extractSDFContours()
    {
        m_contours.clear();
        m_contours.assign(m_unrolledMeshes.size(), zObjectGraph());

        for (int i = 0; i < static_cast<int>(m_unrolledMeshes.size()); ++i)
        {
            zFnMesh fnMesh(m_unrolledMeshes[i]);
            zPoint minBounds;
            zPoint maxBounds;
            fnMesh.getBounds(minBounds, maxBounds);
            z3dpCreateRectangleGraph(m_contours[i], minBounds, maxBounds, std::max(m_settings.printWidth, 0.0f));
        }
    }

    zObjectMeshArray& z3DPSlicer::sliceMeshes() { return m_sliceMeshes; }
    zObjectMeshArray& z3DPSlicer::triangulatedMeshes() { return m_triangulatedMeshes; }
    zObjectMeshArray& z3DPSlicer::unrolledMeshes() { return m_unrolledMeshes; }
    zObjectMeshScalarFieldArray& z3DPSlicer::fields() { return m_fields; }
    zObjectGraphArray& z3DPSlicer::contours() { return m_contours; }
    zObjectGraphArray& z3DPSlicer::edgeLoops() { return m_features.edgeLoops; }
    zObjectPointCloud& z3DPSlicer::cornerPoints() { return m_features.cornerPoints; }
    zObjectMesh& z3DPSlicer::topMesh() { return m_features.topMesh; }
    zObjectMesh& z3DPSlicer::bottomMesh() { return m_features.bottomMesh; }

    const zObjectMeshArray& z3DPSlicer::sliceMeshes() const { return m_sliceMeshes; }
    const zObjectMeshArray& z3DPSlicer::unrolledMeshes() const { return m_unrolledMeshes; }
    const zObjectGraphArray& z3DPSlicer::contours() const { return m_contours; }
    const z3DPSlicingFeatures& z3DPSlicer::features() const { return m_features; }
}
