#include <algorithm>
#include <cmath>

#include <zspace/zToolsets/fabrication/z3DPPrintSynthesis.h>

namespace zSpace
{
    namespace
    {
        zPoint z3dpInterpolate(const zPoint& a, const zPoint& b, float t)
        {
            return zPoint(
                a.x + ((b.x - a.x) * t),
                a.y + ((b.y - a.y) * t),
                a.z + ((b.z - a.z) * t));
        }

        float z3dpDistance(const zPoint& a, const zPoint& b)
        {
            const float dx = b.x - a.x;
            const float dy = b.y - a.y;
            const float dz = b.z - a.z;
            return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
        }

        void z3dpResampleGraph(zObjectGraph& inGraph, zObjectGraph& outGraph, float spacing)
        {
            zFnGraph fnIn(inGraph);
            zPointArray positions;
            zIntArray edges;
            fnIn.getVertexPositions(positions);
            fnIn.getEdgeData(edges);

            zPointArray outPositions;
            zIntArray outEdges;
            const float safeSpacing = std::max(spacing, 0.0001f);

            for (int i = 0; i + 1 < static_cast<int>(edges.size()); i += 2)
            {
                const zPoint& a = positions[edges[i]];
                const zPoint& b = positions[edges[i + 1]];
                const float length = z3dpDistance(a, b);
                const int segments = std::max(1, static_cast<int>(std::ceil(length / safeSpacing)));

                const int start = static_cast<int>(outPositions.size());
                for (int j = 0; j <= segments; ++j)
                {
                    if (i > 0 && j == 0) continue;
                    outPositions.push_back(z3dpInterpolate(a, b, static_cast<float>(j) / static_cast<float>(segments)));
                }

                const int count = static_cast<int>(outPositions.size()) - start;
                for (int j = 0; j < count - 1; ++j)
                {
                    outEdges.push_back(start + j);
                    outEdges.push_back(start + j + 1);
                }
            }

            if (outPositions.size() > 2)
            {
                outEdges.push_back(static_cast<int>(outPositions.size()) - 1);
                outEdges.push_back(0);
            }

            zFnGraph fnOut(outGraph);
            fnOut.clear();
            fnOut.create(outPositions, outEdges);
        }
    }

    void z3DPPrintSynthesis::setSettings(const z3DPSettings& settings)
    {
        m_settings = settings;
    }

    void z3DPPrintSynthesis::setInputs(zObjectGraphArray& contours, zObjectMeshArray& unrolledMeshes, zObjectMeshArray& sliceMeshes)
    {
        m_contours = &contours;
        m_unrolledMeshes = &unrolledMeshes;
        m_sliceMeshes = &sliceMeshes;
    }

    void z3DPPrintSynthesis::resampleContours()
    {
        if (!m_contours) return;
        m_resampledContours.clear();
        m_resampledContours.assign(m_contours->size(), zObjectGraph());

        for (int i = 0; i < static_cast<int>(m_contours->size()); ++i)
        {
            z3dpResampleGraph((*m_contours)[i], m_resampledContours[i], m_settings.printSpacing);
        }
    }

    void z3DPPrintSynthesis::identifyFeaturePoints()
    {
        m_featurePoints.clear();
        m_featurePoints.assign(m_resampledContours.size(), zObjectPointCloud());

        for (int i = 0; i < static_cast<int>(m_resampledContours.size()); ++i)
        {
            zFnGraph fnGraph(m_resampledContours[i]);
            zPointArray positions;
            fnGraph.getVertexPositions(positions);

            zPointArray features;
            if (!positions.empty()) features.push_back(positions.front());
            if (positions.size() > 2) features.push_back(positions[positions.size() / 2]);

            zFnPointCloud fnPoints(m_featurePoints[i]);
            fnPoints.create(features);
        }
    }

    void z3DPPrintSynthesis::mapContoursToSlices()
    {
        m_printPaths = m_resampledContours;

        if (!m_sliceMeshes) return;
        for (int i = 0; i < static_cast<int>(m_printPaths.size()) && i < static_cast<int>(m_sliceMeshes->size()); ++i)
        {
            zFnMesh fnSlice((*m_sliceMeshes)[i]);
            zPoint minBounds;
            zPoint maxBounds;
            fnSlice.getBounds(minBounds, maxBounds);

            zFnGraph fnPath(m_printPaths[i]);
            zPointArray positions;
            fnPath.getVertexPositions(positions);
            for (auto& p : positions) p.z = minBounds.z;
            fnPath.setVertexPositions(positions);
        }
    }

    void z3DPPrintSynthesis::sequenceToolpaths()
    {
    }

    void z3DPPrintSynthesis::buildPrintMesh()
    {
        if (!m_sliceMeshes || m_sliceMeshes->empty()) return;
        m_printMesh = (*m_sliceMeshes)[0];
    }

    zObjectGraphArray& z3DPPrintSynthesis::printPaths() { return m_printPaths; }
    zObjectPointCloudArray& z3DPPrintSynthesis::featurePoints() { return m_featurePoints; }
    zObjectMesh& z3DPPrintSynthesis::printMesh() { return m_printMesh; }

    const zObjectGraphArray& z3DPPrintSynthesis::printPaths() const { return m_printPaths; }
    const zObjectMesh& z3DPPrintSynthesis::printMesh() const { return m_printMesh; }
}

