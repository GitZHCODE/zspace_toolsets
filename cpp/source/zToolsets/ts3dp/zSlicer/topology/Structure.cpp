#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/ts3dp/zSlicer/TopologySlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <deque>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/SparseCholesky>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- PRIVATE HELPERS
//--------------------------

namespace {

// Boundary loop extraction and in-place vertex alignment
void _alignBoundaryVertices(zObjMesh &mesh,
                              std::vector<int> &outOldVertexIds,
                              int startVertexId = 0) {
  outOldVertexIds.clear();
  zFnMesh fn(mesh);
  if (fn.numVertices() == 0 || startVertexId < 0 ||
      startVertexId >= fn.numVertices())
    return;

  // 1. Locate start boundary halfedge
  zItMeshVertex vStart(mesh, startVertexId);
  zItMeshHalfEdgeArray hEdges;
  vStart.getConnectedHalfEdges(hEdges);
  auto it = std::find_if(hEdges.begin(), hEdges.end(),
                         [](zItMeshHalfEdge &he) { return he.onBoundary(); });
  if (it == hEdges.end())
    return;

  zItMeshHalfEdge heStart = *it, heWalk = heStart;
  const int startVID = heStart.getStartVertex().getId();
  zPointArray positions = {heStart.getStartVertex().getPosition()};
  std::vector<int> loopVertexIds = {startVID};

  // 2. Walk unidirectional closed halfedge loop
  for (int step = 0; step < fn.numEdges() + 10; step++) {
    zItMeshVertex endV = heWalk.getVertex();
    if (endV.getId() == startVID)
      break;
    positions.push_back(endV.getPosition());
    loopVertexIds.push_back(endV.getId());
    heWalk = heWalk.getNext();
  }

  const int N = int(positions.size());
  if (N < 3)
    return;

  // 3. Align orientation (counter-clockwise)
  if (signedArea2D(positions) < 0.0) {
    std::reverse(positions.begin() + 1, positions.end());
    std::reverse(loopVertexIds.begin() + 1, loopVertexIds.end());
  }

  outOldVertexIds = loopVertexIds;

  // 4. In-place mesh vertex reordering (boundary vertices at first N indices)
  const int totalVerts = fn.numVertices();
  zPointArray oldPositions;
  fn.getVertexPositions(oldPositions);
  zIntArray polyConnects, polyCounts;
  fn.getPolygonData(polyConnects, polyCounts);

  std::vector<int> oldToNew(totalVerts, -1);
  for (int i = 0; i < N; i++)
    oldToNew[loopVertexIds[i]] = i;
  int nextId = N;
  for (int v = 0; v < totalVerts; v++)
    if (oldToNew[v] < 0)
      oldToNew[v] = nextId++;

  zPointArray alignedPos(totalVerts);
  for (int oldId = 0; oldId < totalVerts; oldId++)
    alignedPos[oldToNew[oldId]] = oldPositions[oldId];
  for (int &id : polyConnects)
    id = oldToNew[id];

  fn.clear();
  fn.create(alignedPos, polyCounts, polyConnects);
}

// The face on edge a-b of the input mesh that lies toward the next V-loop
// (holds one of its vertices); else a Side-tagged one of the two; -1 if a-b
// is not an edge.
inline int _faceTowardColumn(zObjMesh &mesh, int a, int b,
                             const std::unordered_set<int> &nextColumn,
                             const zIntArray &faceTags) {
  zItMeshHalfEdgeArray hes;
  zItMeshVertex(mesh, a).getConnectedHalfEdges(hes);
  for (auto &he : hes) {
    if (he.getVertex().getId() != b)
      continue;
    int candidates[2] = {-1, -1};
    if (!he.onBoundary())
      candidates[0] = he.getFace().getId();
    zItMeshHalfEdge sym = he.getSym();
    if (!sym.onBoundary())
      candidates[1] = sym.getFace().getId();
    for (int face : candidates) {
      if (face < 0)
        continue;
      zIntArray faceVerts;
      zItMeshFace(mesh, face).getVertices(faceVerts);
      for (int v : faceVerts)
        if (nextColumn.count(v))
          return face;
    }
    for (int face : candidates)
      if (face >= 0 && face < int(faceTags.size()) &&
          faceTags[face] == static_cast<int>(VertexTag::Side))
        return face;
    return (candidates[0] >= 0) ? candidates[0] : candidates[1];
  }
  return -1;
}

constexpr int kTagBottom = static_cast<int>(VertexTag::Bottom);
constexpr int kTagTop = static_cast<int>(VertexTag::Top);

// True when he runs inside a floor-cap face.
inline bool _isCapHalfEdge(zItMeshHalfEdge he, const std::vector<bool> &capFace) {
  return !he.onBoundary() && capFace[he.getFace().getId()];
}

// Floor-cap start halfedge, from geometry alone: the cap corner (vertex in
// exactly one cap face) nearest the world origin - or, if the cap has no
// such corner, the nearest cap-boundary vertex - with ties going to the
// lower vertex id; then the cap-boundary edge there leading to another cap
// corner (ties: shorter, then lower end id), oriented to run inside its cap
// face.
bool _findCapStart(zObjMesh &mesh, const std::vector<bool> &capFace,
                   const zIntArray &capFaceCount, zItMeshHalfEdge &outStart) {
  zFnMesh fn(mesh);
  const int numVerts = fn.numVertices();

  const auto isCapBoundaryEdge = [&](zItMeshHalfEdge he) {
    return _isCapHalfEdge(he, capFace) != _isCapHalfEdge(he.getSym(), capFace);
  };
  const auto onCapBoundary = [&](int v) {
    zItMeshHalfEdgeArray hes;
    zItMeshVertex(mesh, v).getConnectedHalfEdges(hes);
    for (auto &he : hes)
      if (isCapBoundaryEdge(he))
        return true;
    return false;
  };

  int start = -1;
  double bestDistSq = (std::numeric_limits<double>::max)();
  for (int pass = 0; pass < 2 && start < 0; pass++) {
    for (int v = 0; v < numVerts; v++) {
      const bool candidate = (pass == 0)
                                 ? capFaceCount[v] == 1
                                 : (capFaceCount[v] > 0 && onCapBoundary(v));
      if (!candidate)
        continue;
      const zPoint pos = zItMeshVertex(mesh, v).getPosition();
      const double distSq = double(pos.x) * pos.x + double(pos.y) * pos.y +
                            double(pos.z) * pos.z;
      if (distSq < bestDistSq) {
        bestDistSq = distSq;
        start = v;
      }
    }
  }
  if (start < 0)
    return false;

  zItMeshHalfEdgeArray hes;
  zItMeshVertex(mesh, start).getConnectedHalfEdges(hes);
  bool found = false, bestCorner = false;
  double bestLength = 0.0;
  int bestEnd = -1;
  zItMeshHalfEdge best;
  for (auto &he : hes) {
    if (!isCapBoundaryEdge(he))
      continue;
    const int end = he.getVertex().getId();
    const bool corner = capFaceCount[end] == 1;
    const double length = he.getLength();
    const bool better =
        !found || (corner && !bestCorner) ||
        (corner == bestCorner &&
         (length < bestLength - 1e-9 ||
          (std::abs(length - bestLength) <= 1e-9 && end < bestEnd)));
    if (better) {
      found = true;
      best = he;
      bestCorner = corner;
      bestLength = length;
      bestEnd = end;
    }
  }
  if (!found)
    return false;
  outStart = _isCapHalfEdge(best, capFace) ? best : best.getSym();
  return true;
}

// The floor cap as its own mesh: faces in breadth-first order from capStart,
// crossing the opposite edge first (for a quad strip this is exactly the
// edge-ring walk), face vertices listed from each face's entry halfedge.
// outOriginalIds: input-mesh id of every cap-mesh vertex. Returns the number
// of cap faces reached.
int _buildCapMesh(zObjMesh &mesh, const std::vector<bool> &capFace,
                  zItMeshHalfEdge capStart, zObjMesh &outCapMesh,
                  zIntArray &outOriginalIds, int &outStartVertex) {
  zFnMesh fn(mesh);
  zPointArray positions;
  zIntArray counts, connects;
  zIntArray vertexMap(fn.numVertices(), -1);
  outOriginalIds.clear();

  std::vector<bool> visited(fn.numPolygons(), false);
  std::deque<zItMeshHalfEdge> queue = {capStart};
  visited[capStart.getFace().getId()] = true;
  int faceCount = 0;
  while (!queue.empty()) {
    zItMeshHalfEdge he = queue.front();
    queue.pop_front();

    zIntArray faceVerts;
    getFaceVerticesFromHalfedge(he, false, faceVerts);
    counts.push_back(int(faceVerts.size()));
    for (int v : faceVerts) {
      int &mapped = vertexMap[v];
      if (mapped < 0) {
        mapped = int(positions.size());
        positions.push_back(zItMeshVertex(mesh, v).getPosition());
        outOriginalIds.push_back(v);
      }
      connects.push_back(mapped);
    }
    faceCount++;

    const zItMeshHalfEdge crossings[4] = {he.getNext().getNext(), he.getNext(),
                                          he.getPrev(), he};
    for (zItMeshHalfEdge crossing : crossings) {
      zItMeshHalfEdge entry = crossing.getSym();
      if (_isCapHalfEdge(entry, capFace) && !visited[entry.getFace().getId()]) {
        visited[entry.getFace().getId()] = true;
        queue.push_back(entry);
      }
    }
  }

  zFnMesh fnCap(outCapMesh);
  fnCap.clear();
  fnCap.create(positions, counts, connects);
  outStartVertex = vertexMap[capStart.getStartVertex().getId()];
  return faceCount;
}

// One column: from floor vertex v along its single edge leaving the floor,
// then straight up the quad edge loop until the first Top vertex.
bool _traceColumn(zObjMesh &mesh, const zIntArray &tags, int v, int maxSteps,
                  zPointArray &outPositions, zScalarArray &outArcLength,
                  zIntArray &outIds, std::string &outError) {
  zItMeshVertex vertex(mesh, v);
  zItMeshHalfEdgeArray hes;
  vertex.getConnectedHalfEdges(hes);
  zItMeshHalfEdge up;
  int upCount = 0;
  for (auto &he : hes) {
    if (tags[he.getVertex().getId()] != kTagBottom) {
      up = he;
      upCount++;
    }
  }
  if (upCount != 1) {
    outError = "floor vertex " + std::to_string(v) + " has " +
               std::to_string(upCount) +
               " edges leaving the floor (needs exactly 1).";
    return false;
  }

  outPositions = {vertex.getPosition()};
  outArcLength = {0.0};
  outIds = {v};
  double accumulatedDist = 0.0;
  zItMeshHalfEdge he = up;
  for (int step = 0; step < maxSteps; step++) {
    zItMeshVertex end = he.getVertex();
    accumulatedDist += he.getLength();
    outPositions.push_back(end.getPosition());
    outArcLength.push_back(accumulatedDist);
    outIds.push_back(end.getId());

    const int tag = tags[end.getId()];
    if (tag == kTagTop)
      return true;
    if (tag == kTagBottom) {
      outError = "the column from floor vertex " + std::to_string(v) +
                 " returns to the floor at vertex " +
                 std::to_string(end.getId()) + ".";
      return false;
    }
    if (end.getValence() != 4) {
      outError = "the column from floor vertex " + std::to_string(v) +
                 " hits irregular vertex " + std::to_string(end.getId()) +
                 " (valence " + std::to_string(end.getValence()) +
                 ") before reaching a Top vertex.";
      return false;
    }
    he = nextEdgeLoop(he);
  }
  outError = "the column from floor vertex " + std::to_string(v) +
             " never reaches a Top vertex.";
  return false;
}

// Uniform-weight harmonic interpolation system for the cap mesh's interior
// vertices (indices >= numBoundary): outL * dInterior = outC * dBoundary.
void _buildInteriorSystem(zObjMesh &capMesh, int numBoundary,
                          Eigen::SparseMatrix<double> &outL,
                          Eigen::SparseMatrix<double> &outC) {
  zFnMesh fn(capMesh);
  const int numInterior = fn.numVertices() - numBoundary;
  outL.resize(0, 0);
  outC.resize(0, 0);
  if (numInterior <= 0)
    return;

  std::vector<Eigen::Triplet<double>> lTriplets, cTriplets;
  std::vector<double> degree(numInterior, 0.0);
  for (zItMeshHalfEdge he(capMesh); !he.end(); he++) {
    const int a = he.getStartVertex().getId(), b = he.getVertex().getId();
    if (a < numBoundary)
      continue;
    degree[a - numBoundary] += 1.0;
    if (b >= numBoundary)
      lTriplets.emplace_back(a - numBoundary, b - numBoundary, -1.0);
    else
      cTriplets.emplace_back(a - numBoundary, b, 1.0);
  }
  for (int k = 0; k < numInterior; k++)
    lTriplets.emplace_back(k, k, degree[k]);

  outL.resize(numInterior, numInterior);
  outL.setFromTriplets(lTriplets.begin(), lTriplets.end());
  outC.resize(numInterior, numBoundary);
  outC.setFromTriplets(cTriplets.begin(), cTriplets.end());
}

} // namespace

//--------------------------
//---- PUBLIC METHODS: STEP 1 (3D SLICING)
//--------------------------

// Floor cap (U-mesh), U-loop, V-loops with per-vertex values, interior
// interpolation system and bracing parameterization - all from the input
// mesh and its vertex tags (see the class comment).
ZSPACE_INLINE bool TopologySlicer::computeBasicStructure(bool useParallel) {
  const bool success = computeBasicStructureImpl(useParallel);
  if (success) zTsUtils::printSuccess("[TopologySlicer::computeBasicStructure] Success");
  else zTsUtils::printError("[TopologySlicer::computeBasicStructure] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computeBasicStructureImpl(bool useParallel) {
  const auto clearStructure = [&]() {
    _VLoops.clear();
    _VScalar.clear();
    _VLoopSegmentWidthIds.clear();
    _InteriorLaplacian.resize(0, 0);
    _InteriorCoupling.resize(0, 0);
    _BracingParams.clear();
    zFnMesh(_UMesh).clear();
  };
  const auto fail = [&](const std::string &message) {
    clearStructure();
    _Status = "failed: " + message;
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  };

  clearStructure();
  o_IsoMeshBoundary.clear();
  o_IsoMeshBoundary2D.clear();
  o_BracingGraphs2D.clear();
  o_IsoMeshes.clear();
  o_IsoMeshes2D.clear();
  o_PathFields.clear();
  o_WidthFields.clear();
  o_TrackFields.clear();
  o_VertexLayerHeights.clear();

  if (!i_SlotPoint) {
    _Status = "failed: i_SlotPoint is required. Call setSlotAnchor() before computing.";
    printDebug("[computeBasicStructureImpl] " + _Status);
    return false;
  }

  if (!i_BlockMesh)
    return fail("i_BlockMesh is null.");

  zObjMesh &mesh = *i_BlockMesh;
  zFnMesh fn(mesh);
  const int numVerts = fn.numVertices();
  if (numVerts == 0 || fn.numPolygons() == 0)
    return fail("mesh has no vertices or faces.");
  std::string tagError;
  if (!zTsUtils::vertexTagsFromFaceTags(mesh, i_FaceTags, _VertexTags,
                                        &tagError))
    return fail("face tags: " + tagError +
                " - pass one tag per face to setInputMesh(mesh, faceTags).");

  // 1. All-quad check and floor-cap faces (every vertex Bottom).
  zIntArray polyConnects, polyCounts;
  fn.getPolygonData(polyConnects, polyCounts);
  std::vector<bool> capFace(polyCounts.size(), false);
  zIntArray capFaceCount(numVerts, 0);
  int capFaces = 0;
  for (size_t f = 0, offset = 0; f < polyCounts.size();
       offset += polyCounts[f], f++) {
    if (polyCounts[f] != 4)
      return fail("face " + std::to_string(f) + " has " +
                  std::to_string(polyCounts[f]) +
                  " vertices - the topology slicer needs an all-quad mesh.");
    bool allBottom = true;
    for (int k = 0; k < 4; k++)
      allBottom = allBottom &&
                  _VertexTags[polyConnects[offset + k]] == kTagBottom;
    if (allBottom) {
      capFace[f] = true;
      capFaces++;
      for (int k = 0; k < 4; k++)
        capFaceCount[polyConnects[offset + k]]++;
    }
  }
  if (capFaces == 0)
    return fail("no face has all-Bottom vertices - tag the print floor "
                "Bottom.");
  if (std::find(_VertexTags.begin(), _VertexTags.end(), kTagTop) ==
      _VertexTags.end())
    return fail("no Top-tagged vertex - tag the print ceiling Top.");

  // 2. U-mesh: the floor cap, from its geometric start corner.
  zItMeshHalfEdge capStart;
  if (!_findCapStart(mesh, capFace, capFaceCount, capStart))
    return fail("could not find a start corner on the floor cap.");
  zIntArray capOriginalIds;
  int capStartVertex = -1;
  const int reached = _buildCapMesh(mesh, capFace, capStart, _UMesh,
                                    capOriginalIds, capStartVertex);
  if (reached != capFaces)
    return fail("the Bottom faces are not one connected patch (" +
                std::to_string(reached) + " of " + std::to_string(capFaces) +
                " reached from the start corner).");

  // 3. Re-order _UMesh: boundary vertices first, in loop order.
  std::vector<int> loopCapIds;
  _alignBoundaryVertices(_UMesh, loopCapIds, capStartVertex);
  const int numLoopVerts = int(loopCapIds.size());
  if (numLoopVerts < 3)
    return fail("the floor cap boundary has fewer than 3 vertices.");
  int boundaryHalfEdges = 0;
  for (zItMeshHalfEdge he(_UMesh); !he.end(); he++)
    if (he.onBoundary())
      boundaryHalfEdges++;
  if (boundaryHalfEdges != numLoopVerts)
    return fail("the floor cap has more than one boundary loop - only a "
                "single loop is supported.");

  std::vector<zIntArray> vLoopVertexIds;
  zIntArray uMeshVertexIds;

  // Input-mesh id of every re-ordered _UMesh vertex (same order as
  // _alignBoundaryVertices(): loop vertices, then the rest ascending).
  {
    std::vector<bool> onLoop(capOriginalIds.size(), false);
    for (int capId : loopCapIds) {
      uMeshVertexIds.push_back(capOriginalIds[capId]);
      onLoop[capId] = true;
    }
    for (size_t capId = 0; capId < capOriginalIds.size(); capId++)
      if (!onLoop[capId])
        uMeshVertexIds.push_back(capOriginalIds[capId]);
  }

  // 4. V-loops: one column per U-loop vertex, up to the first Top vertex.
  _VLoops.reserve(numLoopVerts);
  _VScalar.reserve(numLoopVerts);
  vLoopVertexIds.reserve(numLoopVerts);
  for (int i = 0; i < numLoopVerts; i++) {
    zPointArray positions;
    zScalarArray arcLength;
    zIntArray ids;
    std::string error;
    if (!_traceColumn(mesh, _VertexTags, uMeshVertexIds[i],
                      fn.numEdges() + 1, positions, arcLength, ids, error))
      return fail(error);

    zIntArray edges;
    edges.reserve(2 * (positions.size() - 1));
    for (int k = 0; k + 1 < int(positions.size()); k++) {
      edges.push_back(k);
      edges.push_back(k + 1);
    }
    zObjGraph vGraph;
    zFnGraph(vGraph).create(positions, edges);
    _VLoops.push_back(vGraph);
    _VScalar.push_back(arcLength);
    vLoopVertexIds.push_back(ids);
  }

  // 4b. Width group of the wall face between V-loop i and i + 1 along every
  //     segment of V-loop i (the boundary is split where it changes).
  _VLoopSegmentWidthIds.assign(numLoopVerts, zIntArray());
  for (int i = 0; i < numLoopVerts; i++) {
    const zIntArray &ids = vLoopVertexIds[i];
    const zIntArray &next = vLoopVertexIds[(i + 1) % numLoopVerts];
    const std::unordered_set<int> nextColumn(next.begin(), next.end());
    zIntArray &groups = _VLoopSegmentWidthIds[i];
    groups.assign(ids.size() > 1 ? ids.size() - 1 : 0, 0);
    for (size_t k = 0; k + 1 < ids.size(); k++) {
      const int face =
          _faceTowardColumn(mesh, ids[k], ids[k + 1], nextColumn, i_FaceTags);
      groups[k] = (face >= 0 && face < int(i_FaceWidthIds.size()))
                      ? i_FaceWidthIds[face]
                      : 0;
    }
  }

  // 6. Interior floor-cap vertices follow the boundary harmonically.
  _buildInteriorSystem(_UMesh, numLoopVerts, _InteriorLaplacian,
                       _InteriorCoupling);

  // 7. Bracing graph parameterized on the floor cap.
  if (i_BracingGraph) {
    parameterizeGraph(*i_BracingGraph, _UMesh, _BracingParams);
  }

  // Validate the columns used for slicing.
  for (const auto &arc : _VScalar) {
    if (arc.empty() || !std::isfinite(arc.back()) || arc.back() <= 1e-3)
      return fail("V-loop length is zero or invalid.");
  }
  if (_VScalar.empty())
    return fail("no V-loops are available for slicing.");
  std::ostringstream out;
  out << "Computed structure: floor cap faces=" << capFaces
      << " U-loop vertices=" << numLoopVerts
      << " interior vertices=" << _InteriorLaplacian.rows()
      << " width groups="
      << zTsUtils::widthGroupCount(i_FaceTags, i_FaceWidthIds);
  _Status = out.str();

  return true;
}

// Interpolates layer section meshes/boundaries/bracing along the V-loops.
ZSPACE_INLINE bool TopologySlicer::computeSlice(bool useParallel) {
  const bool success = computeSliceImpl(useParallel);
  if (success) zTsUtils::printSuccess("[TopologySlicer::computeSlice] Success");
  else zTsUtils::printError("[TopologySlicer::computeSlice] Failed");
  return success;
}

ZSPACE_INLINE bool TopologySlicer::computeSliceImpl(bool useParallel) {
  zObjMeshArray &outSectionMeshes = o_IsoMeshes;
  outSectionMeshes.clear();
  o_IsoMeshBoundary.clear();
  o_VertexLayerHeights.clear();

  if (_VLoops.empty() || _VScalar.empty()) {
    _Status = "failed: No V-loops available.";
    printDebug("[computeSliceImpl] " + _Status);
    return false;
  }
  if (i_Settings.LayerSpacing <= 0.0f) {
    _Status = "failed: LayerSpacing must be > 0.";
    printDebug("[computeSliceImpl] " + _Status);
    return false;
  }

  double shortestLoopLength = (std::numeric_limits<double>::max)();
  for (const auto &arc : _VScalar) {
    if (arc.empty() || !std::isfinite(arc.back()) || arc.back() <= 1e-3) {
      _Status = "failed: V-loop length is zero or invalid.";
      printDebug("[computeSliceImpl] " + _Status);
      return false;
    }
    shortestLoopLength = (std::min)(shortestLoopLength, double(arc.back()));
  }

  // ceil(L / LayerSpacing) + 1 layers (the same rule as HarmonicSlicer, on
  // its own L); each V-loop is divided evenly into them.
  const int totalLayers =
      zTsUtils::layersFromSpacing(shortestLoopLength, i_Settings.LayerSpacing);
  if (totalLayers < 2) {
    _Status = "failed: no valid shortest V-loop length - run computeBasicStructure() first.";
    printDebug("[computeSliceImpl] " + _Status);
    return false;
  }

  // Interior floor-cap vertices (no column) are displaced harmonically with
  // the boundary - factorized once here, solved per layer below.
  const int numBoundary = int(_VLoops.size());
  const int numInterior = int(_InteriorLaplacian.rows());
  zPointArray restPositions;
  zFnMesh(_UMesh).getVertexPositions(restPositions);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> interiorSolver;
  if (numInterior > 0) {
    interiorSolver.compute(_InteriorLaplacian);
    if (interiorSolver.info() != Eigen::Success) {
      _Status = "failed: could not factorize the floor-cap interior system.";
      printDebug("[computeSliceImpl] " + _Status);
      return false;
    }
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

  // 1. Interpolate 3D section mesh and extract boundaries
  // in parallel
  outSectionMeshes.assign(totalLayers, _UMesh);
  o_IsoMeshBoundary.assign(totalLayers, zObjGraph());

#pragma omp parallel for schedule(dynamic) if(useParallel)
  for (int l = 0; l < totalLayers; l++) {
    const double ratio = double(l) / (totalLayers - 1);
    zFnMesh fnMesh(outSectionMeshes[l]);
    zPoint *points = fnMesh.getRawVertexPositions();
    if (points) {
      for (size_t i = 0; i < _VLoops.size() && i < _VScalar.size() &&
                         i < fnMesh.numVertices();
           i++) {
        samplePolyline(_VLoops[i], _VScalar[i], ratio, points[i]);
      }
      if (numInterior > 0) {
        Eigen::MatrixXd boundaryDelta(numBoundary, 3);
        for (int b = 0; b < numBoundary; b++) {
          boundaryDelta(b, 0) = points[b].x - restPositions[b].x;
          boundaryDelta(b, 1) = points[b].y - restPositions[b].y;
          boundaryDelta(b, 2) = points[b].z - restPositions[b].z;
        }
        const Eigen::MatrixXd interiorDelta =
            interiorSolver.solve(_InteriorCoupling * boundaryDelta);
        for (int k = 0; k < numInterior; k++) {
          const zPoint &rest = restPositions[numBoundary + k];
          points[numBoundary + k] = zPoint(rest.x + interiorDelta(k, 0),
                                           rest.y + interiorDelta(k, 1),
                                           rest.z + interiorDelta(k, 2));
        }
      }
    }

    std::vector<int> secOldVertexIds;
    _alignBoundaryVertices(outSectionMeshes[l], secOldVertexIds, 0);
    const int boundaryCount = int(secOldVertexIds.size());
    if (boundaryCount >= 3) {
      zPointArray boundaryPositions;
      fnMesh.getVertexPositions(boundaryPositions);
      boundaryPositions.resize(boundaryCount);
      zIntArray edgeConnects;
      edgeConnects.reserve(boundaryCount * 2);
      for (int i = 0; i < boundaryCount; ++i) {
        edgeConnects.push_back(i);
        edgeConnects.push_back((i + 1) % boundaryCount);
      }
      zFnGraph(o_IsoMeshBoundary[l]).create(boundaryPositions, edgeConnects);
    }

  }

  // Per-vertex layer spacing.
  o_VertexLayerHeights.assign(_VScalar.size(), 0.0f);
  if (totalLayers > 1) {
    const float denom = float(totalLayers - 1);
    for (size_t k = 0; k < _VScalar.size(); k++) {
      if (!_VScalar[k].empty())
        o_VertexLayerHeights[k] = float(_VScalar[k].back() / denom);
    }
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Computed slices: V-loops=" << _VLoops.size()
      << " sectionMeshes=" << outSectionMeshes.size() << " in " << std::fixed
      << std::setprecision(2) << wallTotalMs << " ms";
  _Status = out.str();

  return true;
}

} // namespace zSpace
