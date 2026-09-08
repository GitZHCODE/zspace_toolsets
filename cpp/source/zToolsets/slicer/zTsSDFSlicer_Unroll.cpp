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

#include <headers/zToolsets/slicer/zTsSDFSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

namespace {

// Dual Graph Computation and 2D Unrolling
bool _computeDualGraph(zObjMesh &mesh, zIntPairArray &bsfHalfEdgePairs) {
  bsfHalfEdgePairs.clear();
  zFnMesh fnMesh(mesh);
  if (fnMesh.numPolygons() <= 1)
    return false;

  zObjGraph graph;
  zIntArray inEdgeDualEdge;
  zIntArray dualEdgeInEdge;
  fnMesh.getDualGraph(graph, inEdgeDualEdge, dualEdgeInEdge, true, false,
                      false);

  int maxValence = -1;
  zItGraphVertex maxVertex;
  for (zItGraphVertex v(graph); !v.end(); v++) {
    if (v.getValence() > maxValence) {
      maxValence = v.getValence();
      maxVertex = v;
    }
  }
  if (maxValence < 0)
    return false;

  zItGraphVertexArray bsfVertices;
  zIntPairArray bsfFacePairs;
  maxVertex.getBSF(bsfVertices, bsfFacePairs);
  if (static_cast<int>(bsfFacePairs.size()) + 1 != fnMesh.numPolygons())
    return false;

  bsfHalfEdgePairs.reserve(bsfFacePairs.size());
  for (const auto &facePair : bsfFacePairs) {
    zItMeshFace f1(mesh, facePair.first);
    zItMeshFace f2(mesh, facePair.second);
    zItMeshHalfEdgeArray f1HE, f2HE;
    f1.getHalfEdges(f1HE);
    f2.getHalfEdges(f2HE);

    zIntPair commonHE(-1, -1);
    for (auto &he1 : f1HE) {
      for (auto &he2 : f2HE) {
        if (he1.getEdge().getId() == he2.getEdge().getId()) {
          commonHE = zIntPair(he1.getId(), he2.getId());
          break;
        }
      }
      if (commonHE.first != -1)
        break;
    }

    if (commonHE.first == -1)
      return false;
    bsfHalfEdgePairs.push_back(commonHE);
  }

  return true;
}

bool _unrollMesh(zObjMesh &sourceMesh, const zIntPairArray &bsfHalfEdgePairs,
                 zObjMesh &flattenedMesh) {
  zFnMesh fnSource(sourceMesh);
  const int numVerts = fnSource.numVertices();
  const int numFaces = fnSource.numPolygons();
  if (numVerts == 0 || numFaces == 0)
    return false;

  zPoint *srcPositions = fnSource.getRawVertexPositions();
  if (!srcPositions)
    return false;

  zPointArray flatPositions(numVerts, zPoint(0, 0, 0));
  std::vector<bool> isPlaced(numVerts, false);

  auto placeFirstEdge = [&](zItMeshHalfEdge he) -> bool {
    const float edgeLength = he.getLength();
    if (edgeLength <= 1e-6f)
      return false;
    const int pStartId = he.getStartVertex().getId();
    const int pEndId = he.getVertex().getId();
    flatPositions[pStartId] = zPoint(0, 0, 0);
    flatPositions[pEndId] = zPoint(0, edgeLength, 0);
    isPlaced[pStartId] = true;
    isPlaced[pEndId] = true;
    return true;
  };

  auto unrollFace = [&](zItMeshFace face, zItMeshHalfEdge baseHE) -> bool {
    const float edgeLength = baseHE.getLength();
    const int pStartId = baseHE.getStartVertex().getId();
    const int pEndId = baseHE.getVertex().getId();
    if (edgeLength <= 1e-6f || !isPlaced[pStartId] || !isPlaced[pEndId])
      return false;

    const zPoint A = srcPositions[pStartId], B = srcPositions[pEndId];
    const zPoint a = flatPositions[pStartId], b = flatPositions[pEndId];
    const float denom = edgeLength * edgeLength;
    zVector ba(B.x - A.x, B.y - A.y, B.z - A.z);

    zItMeshHalfEdge walker = baseHE;
    for (int v = 0; v < face.getNumVertices(); v++) {
      walker = walker.getNext();
      const int vid = walker.getVertex().getId();
      if (isPlaced[vid])
        continue;

      const zPoint C = srcPositions[vid];
      zVector ca(C.x - A.x, C.y - A.y, C.z - A.z);
      const float s = static_cast<float>((ba ^ ca).length() / denom);
      const float c = static_cast<float>((ba * ca) / denom);
      flatPositions[vid] = zPoint(a.x + c * (b.x - a.x) - s * (b.y - a.y),
                                  a.y + c * (b.y - a.y) + s * (b.x - a.x), 0.0);
      isPlaced[vid] = true;
    }
    return true;
  };

  if (numFaces == 1) {
    zItMeshFace face(sourceMesh, 0);
    zItMeshHalfEdge he = face.getHalfEdge();
    if (!placeFirstEdge(he) || !unrollFace(face, he))
      return false;
  } else {
    if (bsfHalfEdgePairs.empty())
      return false;

    for (size_t pairId = 0; pairId < bsfHalfEdgePairs.size(); pairId++) {
      zItMeshHalfEdge parentHE(sourceMesh, bsfHalfEdgePairs[pairId].first);
      zItMeshHalfEdge childHE(sourceMesh, bsfHalfEdgePairs[pairId].second);

      if (pairId == 0) {
        if (!placeFirstEdge(parentHE) ||
            !unrollFace(parentHE.getFace(), parentHE))
          return false;
      }

      if (!unrollFace(childHE.getFace(), childHE))
        return false;
    }
  }

  flattenedMesh = sourceMesh;
  zFnMesh fnFlattened(flattenedMesh);
  fnFlattened.setVertexPositions(flatPositions);
  return true;
}
} // anonymous namespace

//---- MAIN COMPUTE METHODS: STEP 2 (2D UNROLLING)

ZSPACE_INLINE bool zTsSDFSlicer::computeUnroll() {
  _2DSectionMeshes.clear();
  _2DSectionGraphs.clear();
  _2DBracingGraphs.clear();

  const int totalLayers = int(sectionMeshes.size());
  if (totalLayers == 0) {
    m_Status = "failed: No section meshes available for unroll.";
    printError("[computeUnroll] " + m_Status);
    return false;
  }

  zIntPairArray bsfHalfEdgePairs;
  if (zFnMesh(sectionMeshes[0]).numPolygons() > 1) {
    if (!_computeDualGraph(sectionMeshes[0], bsfHalfEdgePairs)) {
      m_Status = "failed: Could not compute dual graph for unroll.";
      printError("[computeUnroll] " + m_Status);
      return false;
    }
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

  _2DSectionMeshes.assign(totalLayers, zObjMesh());
  _2DSectionGraphs.assign(totalLayers, zObjGraph());
  if (!m_BracingParams.empty()) {
    _2DBracingGraphs.assign(totalLayers, zObjGraph());
  }

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    if (_unrollMesh(sectionMeshes[l], bsfHalfEdgePairs,
                    _2DSectionMeshes[l])) {
      zFnMesh fnFlatMesh(_2DSectionMeshes[l]);
      zPointArray flatMeshPos;
      fnFlatMesh.getVertexPositions(flatMeshPos);

      zFnGraph fn3DGraph(sectionGraphs[l]);
      const int N = fn3DGraph.numVertices();
      zIntArray edgeConnects;
      fn3DGraph.getEdgeData(edgeConnects);

      if (N > 0 && N <= (int)flatMeshPos.size()) {
        zPointArray flatBoundPositions(flatMeshPos.begin(),
                                       flatMeshPos.begin() + N);
        zFnGraph fnFlatGraph(_2DSectionGraphs[l]);
        fnFlatGraph.create(flatBoundPositions, edgeConnects);

        if (!m_BracingParams.empty()) {
          realizeParameterGraph(
              m_BracingParams, _2DSectionMeshes[l], _2DBracingGraphs[l]);
        }
      }
    }
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Computed unroll: " << _2DSectionMeshes.size() << " layers in "
      << std::fixed << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  printSuccess("[computeUnroll] " + m_Status);
  return true;
}

//---- MAIN COMPUTE METHODS: STEP 5 (PROJECT 2D PRINT PATHS BACK TO 3D)
//---- Paired with computeUnroll() above - unroll flattens 3D -> 2D, this
//---- is the inverse mapping, 2D -> 3D, for the (possibly post-processed)
//---- print paths.

ZSPACE_INLINE bool zTsSDFSlicer::projectPathsTo3D(
    const std::vector<PrintPath> &in2DPaths,
    std::vector<PrintPath> &out3DPaths) {
  const int totalLayers = int(in2DPaths.size());
  if (totalLayers == 0 || _2DSectionMeshes.empty() || sectionMeshes.empty()) {
    m_Status = "failed: No 2D print paths or section meshes available for "
               "projection.";
    printError("[projectPathsTo3D] " + m_Status);
    return false;
  }

  out3DPaths.assign(totalLayers, PrintPath());
  _3DPrintPaths.assign(totalLayers, PrintPath());

  zFloatArray meshVertexHeights;
  if (!m_VScalar.empty() && totalLayers > 1) {
    meshVertexHeights.resize(m_VScalar.size());
    const float denom = float(totalLayers - 1);
    for (size_t k = 0; k < m_VScalar.size(); k++) {
      meshVertexHeights[k] = (!m_VScalar[k].empty())
                                 ? float(m_VScalar[k].back() / denom)
                                 : i_Settings.layerSpacing;
    }
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    out3DPaths[l].graph = in2DPaths[l].graph;
    out3DPaths[l].widths = in2DPaths[l].widths;

    zFnGraph fnContour(out3DPaths[l].graph);
    if (fnContour.numVertices() > 0 && l < int(_2DSectionMeshes.size()) &&
        l < int(sectionMeshes.size())) {
      remapGraph(
          out3DPaths[l].graph, _2DSectionMeshes[l], sectionMeshes[l],
          &out3DPaths[l].normals, &meshVertexHeights, &out3DPaths[l].heights);
    }

    _3DPrintPaths[l] = out3DPaths[l];
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Projected " << out3DPaths.size() << " 3D PrintPaths in "
      << std::fixed << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  printSuccess("[projectPathsTo3D] " + m_Status);
  return true;
}

} // namespace zSpace
