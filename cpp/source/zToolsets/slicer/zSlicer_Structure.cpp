// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Vishu Bhooshan <vishu.bhooshan@zaha-hadid.com>
//

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/slicer/zSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- PRIVATE HELPERS
//--------------------------

namespace {

// UV Start Halfedge Search
bool _getUVStart(zObjMesh &mesh, int startVID, zItMeshHalfEdge &heUStart,
                 zItMeshHalfEdge &heVStart, int &numSegments) {
  zFnMesh fnMesh(mesh);
  zItMeshVertex vStart(mesh, startVID);
  heVStart = zUpEdge(vStart);
  if (heVStart.getVertex().getId() < 0)
    return false;

  const int vId_v = heVStart.getVertex().getId();
  zItMeshHalfEdgeArray hEdgesStart;
  vStart.getConnectedHalfEdges(hEdgesStart);

  bool foundUStart = false;
  for (auto &he : hEdgesStart) {
    if (he != heVStart && he.getVertex().getValence() == 3) {
      heUStart = he;
      foundUStart = true;
      break;
    }
  }

  if (!foundUStart)
    return false;

  zIntArray faceVerts;
  getFaceVerticesFromHalfedge(heUStart, false, faceVerts);
  if (std::find(faceVerts.begin(), faceVerts.end(), vId_v) != faceVerts.end()) {
    heUStart = heUStart.getSym();
  }

  numSegments = 0;
  zItMeshHalfEdge heWalkUp = heVStart;
  while (numSegments < fnMesh.numPolygons() + 10) {
    numSegments++;
    if (heWalkUp.getVertex().getValence() == 3) {
      break;
    }
    heWalkUp = nextEdgeLoop(heWalkUp);
  }

  return numSegments > 0;
}

// Boundary loop extraction, in-place vertex alignment, and corner split detection
void _createBoundaryEdgeGraph(zObjMesh &mesh, zObjGraph &outGraph,
                              std::vector<int> &outOldVertexIds,
                              int startVertexId = 0, bool topoSplit = true,
                              bool angleSplit = true,
                              float angleThresholdDeg = 90.0f,
                              std::vector<int> *outSplitIndices = nullptr) {
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
  std::vector<bool> isCorner = {heStart.getStartVertex().getValence() == 2};

  // 2. Walk unidirectional closed halfedge loop
  for (int step = 0; step < fn.numEdges() + 10; step++) {
    zItMeshVertex endV = heWalk.getVertex();
    if (endV.getId() == startVID)
      break;
    positions.push_back(endV.getPosition());
    loopVertexIds.push_back(endV.getId());
    isCorner.push_back(endV.getValence() == 2);
    heWalk = heWalk.getNext();
  }

  const int N = int(positions.size());
  if (N < 3)
    return;

  // 3. Align orientation (counter-clockwise)
  double signedArea = 0.0;
  for (int i = 0; i < N; i++) {
    int nxt = (i + 1) % N;
    signedArea += (positions[i].x * positions[nxt].y -
                   positions[nxt].x * positions[i].y);
  }
  if (signedArea < 0.0) {
    std::reverse(positions.begin() + 1, positions.end());
    std::reverse(loopVertexIds.begin() + 1, loopVertexIds.end());
    std::reverse(isCorner.begin() + 1, isCorner.end());
  }

  // 4. Construct 2D/3D closed boundary graph
  zIntArray edgeConnects;
  edgeConnects.reserve(N * 2);
  for (int i = 0; i < N; i++) {
    edgeConnects.push_back(i);
    edgeConnects.push_back((i + 1) % N);
  }
  zFnGraph(outGraph).create(positions, edgeConnects);
  outOldVertexIds = loopVertexIds;

  // 5. In-place mesh vertex reordering (boundary vertices at first N indices)
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

  // 6. Corner topological and geometric sharp angle split detection
  if (outSplitIndices) {
    std::set<int> splitSet;
    const double angleRadThresh =
        angleThresholdDeg * 3.14159265358979323846 / 180.0;
    for (int i = 0; i < N; i++) {
      if (topoSplit && isCorner[i])
        splitSet.insert(i);
      if (angleSplit) {
        zVector v1 = vecSub(positions[(i - 1 + N) % N], positions[i]);
        zVector v2 = vecSub(positions[(i + 1) % N], positions[i]);
        if (v1.length() > 1e-6 && v2.length() > 1e-6) {
          v1.normalize();
          v2.normalize();
          double dot = (std::max)(-1.0, (std::min)(1.0, vecDot(v1, v2)));
          if (std::acos(dot) < angleRadThresh)
            splitSet.insert(i);
        }
      }
    }
    outSplitIndices->assign(splitSet.begin(), splitSet.end());
  }
}

// Extracts U-mesh, U-loop, and V-loop graphs by walking the UV grid.
bool _walkUVLoop(zObjMesh &mesh, zItMeshHalfEdge heUStart,
                 zItMeshHalfEdge heVStart, int numSegments, bool topoSplit,
                 bool angleSplit, float angleThresholdDeg, zObjMesh &uMeshObj,
                 zObjGraph &uLoopObj, zObjGraphArray &vLoopGraphs,
                 std::vector<zScalarArray> &vScalars,
                 std::vector<int> &outSplitIndices) {
  zFnMesh fn(mesh);
  vLoopGraphs.clear();
  vScalars.clear();
  outSplitIndices.clear();
  zFnMesh(uMeshObj).clear();
  zFnGraph(uLoopObj).clear();

  auto stripReachedCorner = [](zItMeshHalfEdge he) {
    return he.getStartVertex().getValence() == 3 ||
           he.getVertex().getValence() == 3;
  };

  zPointArray bottomPositions;
  zIntArray bottomCounts;
  zIntArray bottomConnects;
  zIntArray bottomVertexMap(fn.numVertices(), -1);
  std::vector<int> bottomOriginalVertexIds;

  auto appendFaceToMesh = [&](zItMeshHalfEdge he) {
    zIntArray faceVerts;
    getFaceVerticesFromHalfedge(he, false, faceVerts);
    bottomCounts.push_back(int(faceVerts.size()));

    for (int origVID : faceVerts) {
      int &mappedId = bottomVertexMap[origVID];
      if (mappedId < 0) {
        mappedId = int(bottomPositions.size());
        bottomPositions.push_back(zItMeshVertex(mesh, origVID).getPosition());
        bottomOriginalVertexIds.push_back(origVID);
      }
      bottomConnects.push_back(mappedId);
    }
  };

  zItMeshHalfEdge heUWalk = heUStart;
  for (int safety = 0; safety < fn.numPolygons() + 10; safety++) {
    appendFaceToMesh(heUWalk);

    zItMeshHalfEdge nextU = nextEdgeRing(heUWalk);
    zIntArray nextFaceVerts;
    getFaceVerticesFromHalfedge(nextU, false, nextFaceVerts);
    if (nextFaceVerts.size() != 4 || stripReachedCorner(nextU)) {
      break;
    }

    heUWalk = nextU;
  }

  if (bottomPositions.empty())
    return false;

  zFnMesh fnU(uMeshObj);
  fnU.clear();
  fnU.create(bottomPositions, bottomCounts, bottomConnects);

  const int uMeshStartVID = bottomVertexMap[heUStart.getStartVertex().getId()];
  std::vector<int> oldUMeshVertexIds;
  _createBoundaryEdgeGraph(uMeshObj, uLoopObj, oldUMeshVertexIds, uMeshStartVID,
                           topoSplit, angleSplit, angleThresholdDeg,
                           &outSplitIndices);

  zFnGraph fnULoop(uLoopObj);

  const int numLoopVerts = fnULoop.numVertices();
  vLoopGraphs.reserve(numLoopVerts);
  vScalars.reserve(numLoopVerts);

  for (int i = 0; i < numLoopVerts; i++) {
    zItMeshHalfEdge heUp;

    const int oldVID_uMesh = oldUMeshVertexIds[i];
    const int origVID = bottomOriginalVertexIds[oldVID_uMesh];
    if (origVID < 0)
      continue;

    heUp = zUpEdge(zItMeshVertex(mesh, origVID));
    if (heUp.getVertex().getId() < 0)
      continue;

    zPointArray vPositions;
    zIntArray vEdges;
    zScalarArray scalarLoop;
    vPositions.reserve(numSegments + 1);
    vEdges.reserve(numSegments * 2);
    scalarLoop.reserve(numSegments + 1);

    vPositions.push_back(heUp.getStartVertex().getPosition());
    scalarLoop.push_back(0.0);

    double accumulatedDist = 0.0;
    zItMeshHalfEdge heWalk = heUp;
    for (int s = 0; s < numSegments; s++) {
      const zPoint pNext = heWalk.getVertex().getPosition();
      vPositions.push_back(pNext);
      vEdges.push_back(s);
      vEdges.push_back(s + 1);

      accumulatedDist += heWalk.getLength();
      scalarLoop.push_back(accumulatedDist);

      if (s + 1 < numSegments) {
        heWalk = nextEdgeLoop(heWalk);
      }
    }

    zObjGraph vGraph;
    zFnGraph fnV(vGraph);
    fnV.create(vPositions, vEdges);
    vLoopGraphs.push_back(vGraph);
    vScalars.push_back(scalarLoop);
  }

  return !vLoopGraphs.empty();
}

} // anonymous namespace

//--------------------------
//---- PUBLIC METHODS: STEP 1 (3D SLICING)
//--------------------------

// Extracts corner UV start, U-mesh, U-loop, V-loops, and bracing parameterization.
ZSPACE_INLINE bool zSlicer::computeBasicStructure() {
  m_VLoops.clear();
  m_VScalar.clear();
  sectionMeshes.clear();
  sectionGraphs.clear();
  bracingGraphs.clear();
  sectionMeshes2D.clear();
  sectionGraphs2D.clear();
  bracingGraphs2D.clear();
  pathFields.clear();
  widthFields.clear();
  trimSlotGraphs2D.clear();
  m_BracingParams.clear();
  zFnMesh(m_UMesh).clear();
  zFnGraph(m_ULoop).clear();
  zFnMeshScalarField(m_ScalarField).clear();

  if (!i_InputMesh) {
    m_Status = "failed: i_InputMesh is null.";
    printError("[computeBasicStructure] " + m_Status);
    return false;
  }

  zFnMesh fn(*i_InputMesh);
  if (fn.numVertices() == 0 || fn.numPolygons() == 0) {
    m_Status = "failed: mesh has no vertices or faces.";
    printError("[computeBasicStructure] " + m_Status);
    return false;
  }

  const int startVID = (i_Settings.StartCornerVertexID >= 0 &&
                        i_Settings.StartCornerVertexID < fn.numVertices())
                           ? i_Settings.StartCornerVertexID
                           : 0;

  zItMeshHalfEdge heUStart, heVStart;
  int numSegments = 0;
  if (!_getUVStart(*i_InputMesh, startVID, heUStart, heVStart, numSegments)) {
    m_Status = "failed: could not find valid UV start halfedge.";
    printError("[computeBasicStructure] " + m_Status);
    return false;
  }

  m_BoundarySplitIndices.clear();

  if (!_walkUVLoop(*i_InputMesh, heUStart, heVStart, numSegments,
                   i_Settings.TopoSplit, i_Settings.AngleSplit,
                   i_Settings.AngleSplitTol, m_UMesh, m_ULoop,
                   m_VLoops, m_VScalar, m_BoundarySplitIndices)) {
    m_VLoops.clear();
    m_VScalar.clear();
    m_BoundarySplitIndices.clear();
    zFnMesh(m_UMesh).clear();
    zFnGraph(m_ULoop).clear();
    m_Status = "failed: could not extract Umesh and V-loops.";
    printError("[computeBasicStructure] " + m_Status);
    return false;
  }

  // Compute barycentric parameterization of bracingGraph on 3D base mesh
  if (i_BracingGraph) {
    parameterizeGraph(*i_BracingGraph, m_UMesh, m_BracingParams);
  }

  return true;
}

// Interpolates layer section meshes/boundaries/bracing along the V-loops.
ZSPACE_INLINE bool zSlicer::computeSlice() {
  sectionMeshes.clear();
  sectionGraphs.clear();
  bracingGraphs.clear();

  if (m_VLoops.empty() || m_VScalar.empty() ||
      i_Settings.LayerSpacing <= 0.0f) {
    m_Status = "failed: No V-loops available.";
    printError("[computeSlice] " + m_Status);
    return false;
  }

  float minLoopLength = (std::numeric_limits<float>::max)();
  for (const auto &sLoop : m_VScalar) {
    if (!sLoop.empty() && sLoop.back() > 1e-6) {
      minLoopLength = (std::min)(minLoopLength, float(sLoop.back()));
    }
  }

  if (minLoopLength == (std::numeric_limits<float>::max)() ||
      minLoopLength <= 1e-6f) {
    sectionMeshes.clear();
    m_Status = "failed: V-loops length is zero or invalid.";
    printError("[computeSlice] " + m_Status);
    return false;
  }

  const int totalContours =
      (std::max)(1, int(std::ceil(minLoopLength / i_Settings.LayerSpacing)));
  const int totalLayers = totalContours + 1;

  const auto wallStart = std::chrono::high_resolution_clock::now();

  // 1. Interpolate 3D section mesh, extract boundaries, and reconstruct bracing in parallel
  sectionMeshes.assign(totalLayers, m_UMesh);
  sectionGraphs.assign(totalLayers, zObjGraph());
  if (!m_BracingParams.empty()) {
    bracingGraphs.assign(totalLayers, zObjGraph());
  }

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    const double ratio = double(l) / (totalLayers - 1);
    zFnMesh fnMesh(sectionMeshes[l]);
    zPoint *points = fnMesh.getRawVertexPositions();
    if (points) {
      for (size_t i = 0; i < m_VLoops.size() && i < m_VScalar.size() &&
                         i < fnMesh.numVertices();
           i++) {
        samplePolyline(m_VLoops[i], m_VScalar[i], ratio, points[i]);
      }
    }

    std::vector<int> secOldVertexIds;
    _createBoundaryEdgeGraph(sectionMeshes[l], sectionGraphs[l],
                             secOldVertexIds, 0);

    if (!m_BracingParams.empty()) {
      realizeParameterGraph(
          m_BracingParams, sectionMeshes[l], bracingGraphs[l]);
    }
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Computed slices: V-loops=" << m_VLoops.size()
      << " sectionMeshes=" << sectionMeshes.size() << " in " << std::fixed
      << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  printSuccess("[computeSlice] " + m_Status);
  return true;
}

} // namespace zSpace
