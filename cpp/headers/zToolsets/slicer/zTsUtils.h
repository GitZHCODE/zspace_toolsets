// This file is part of zspace, a simple C++ collection of geometry
// data-structures & algorithms, data analysis & visualization framework.
//
// Copyright (C) 2019 ZSPACE
//
// This Source Code Form is subject to the terms of the MIT License
// If a copy of the MIT License was not distributed with this file, You can
// obtain one at https://opensource.org/licenses/MIT.
//
// Author : Weijie Tao <weijie.tao@zspace.com>
//

#ifndef ZSPACE_TS_SLICER_UTILS_H
#define ZSPACE_TS_SLICER_UTILS_H

#pragma once

#include <headers/base/zSpace_Toolsets.h>

#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>

#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Collection of toolsets for geometry, analysis, and fabrication.
 *  @{
 */

/** \addtogroup zTsGeometry
 *	\brief Toolsets for geometry generation, slicing, and manipulation.
 *  @{
 */

//--------------------------
//---- COMMON CONTAINER TYPEDEFS
//--------------------------

typedef std::vector<zObjMesh> zObjMeshArray;
typedef std::vector<zObjGraph> zObjGraphArray;
typedef std::vector<zObjGraph *> zObjGraphPointerArray;
// NOTE: zObjMeshScalarFieldArray is intentionally not declared here. This
// header does not include zObjMeshField.h (nothing in the shared-utils layer
// needs zObjMeshScalarField), so redeclaring the typedef against an
// undeclared type here fails to compile, and once zObjMeshField.h's own
// `using zObjMeshScalarFieldArray = ...` is later seen (via zTsSDFSlicer.h),
// MSVC reports it as a redefinition with a mismatched type. Consumers that
// need the array type should include zObjMeshField.h directly.

//--------------------------
//---- COMMON DATA STRUCTURES
//--------------------------

/*! \struct SlicerSettings
 *	\brief Configuration settings for geodesic slicing and SDF
 * field/toolpath generation. Consumed only by zTsSDFSlicer.
 */
struct SlicerSettings {
  int cornerStartVertexId = 21;
  float layerSpacing = 1.0f;
  bool topoSplit = true;
  bool angleSplit = true;
  float geometricAngleThreshold = 90.0f; // degrees (corner angle < threshold)

  // SDF & Toolpath Parameters
  bool doublelayerboundary = false;
  std::vector<float> printBoundaryWidth = {0.5f};
  float printBracingWidth = 0.5f;
  float printOverlapWidth = 0.0f;
  float trimSlotStaggerEven = 0.5f;
  float trimSlotStaggerOdd = 0.5f;
  int boundaryTrimSegmentId = 1;
  float boundaryTrimSlotRatioEven = 0.45f;
  float boundaryTrimSlotRatioOdd = 0.55f;

  float fieldBoundingScale = 1.05f;
  float fieldCellSize = 0.05f; // Target square grid cell size (default 0.05)
};

/*! \struct PostProcessSettings
 *	\brief Configuration settings for 2D toolpath cleanup and resampling.
 * Consumed only by zTsPrintSyn.
 */
struct PostProcessSettings {
  float sampleLength = 0.15f; // arc-length resampling step (mm)
  float featureAngleThreshold =
      90.0f; // geometric feature corner angle threshold (deg)
};

/*! \struct BracingGraph
 *	\brief Universal graph representation parameterized by 2D/3D mesh face
 * barycentric coordinates.
 */
struct BracingGraph {
  struct Vertex {
    int v0, v1, v2; // Triangle 3 vertex indices
    float u, v, w;  // Barycentric coordinate weights
  };

  std::vector<Vertex> vertices; // Parameterized vertices array
  zIntArray edgeConnects; // Edge topology connectivity [e0_v0, e0_v1, e1_v0,
                          // e1_v1, ...]

  bool empty() const { return vertices.empty() || edgeConnects.empty(); }
  void clear() {
    vertices.clear();
    edgeConnects.clear();
  }
};

/*! \struct PrintPath
 *	\brief Encapsulates a complete printing toolpath including its graph,
 * line widths, layer heights, and normals.
 *	\details Shared exchange type between zTsSDFSlicer (produces raw 2D/3D
 * paths) and zTsPrintSyn (cleans up and resamples paths); neither engine class
 * depends on the other.
 */
struct PrintPath {
  zObjGraph graph;
  zFloatArray widths;
  zFloatArray heights;
  zVectorArray normals;
};

//--------------------------
//---- COMMON GEOMETRY / TOPOLOGY UTILITIES
//--------------------------
// Shared low-level building blocks used by both zTsSDFSlicer and zTsPrintSyn.
// Declared here, defined in zTsUtils.cpp (see ZSPACE_TOOLSETS_INLINE below).

namespace zTsUtils {

//---- Basic Vector Utilities

zVector vecSub(const zPoint &a, const zPoint &b);
double vecDistSq(const zPoint &a, const zPoint &b);
double vecDot(const zVector &a, const zVector &b);
zVector vecCross(const zVector &a, const zVector &b);
unsigned long long makeEdgeKey(int u, int v);

//---- HalfEdge Mesh Topology Walk Utilities

zItMeshHalfEdge nextEdgeRing(zItMeshHalfEdge he);
zItMeshHalfEdge nextEdgeLoop(zItMeshHalfEdge he);
zItMeshHalfEdge zUpEdge(zItMeshVertex v);
void getFaceVerticesFromHalfedge(zItMeshHalfEdge &heStart, bool forward,
                                 zIntArray &faceVerts);

//---- Combine Multiple Graph Objects

void combineGraphObjects(const zObjGraphArray &sources, zObjGraph &outGraph);

//---- Universal Polyline Sampling by Arc-Length with Optional Tangent Vector

bool samplePolyline(const zPoint *positions, const zScalar *scalars,
                    size_t count, double t, zPoint &outPoint,
                    zVector *outTangent = nullptr);
bool samplePolyline(const zObjGraph &polylineGraph, const zScalarArray &scalars,
                    double t, zPoint &outPoint, zVector *outTangent = nullptr);

//---- Universal Mesh Barycentric Parameterizer (2D flat and 3D space meshes)

bool parameterizeGraphOnMesh(const zObjGraph &inGraph, zObjMesh &inMesh,
                             BracingGraph &outParams);

//---- Universal Barycentric Mesh Graph Reconstructor

bool reconstructGraphFromBarycentricParams(
    const BracingGraph &inParams, zObjMesh &targetMesh,
    zObjGraph &outGraph, zVectorArray *outNormals = nullptr,
    const zFloatArray *vertexScalars = nullptr,
    zFloatArray *outHeights = nullptr);

//---- Barycentric Mesh-to-Mesh Graph Projection

bool barycentricProjectionTriMesh(zObjGraph &graph, zObjMesh &inMesh2D,
                                  zObjMesh &projectionMesh3D,
                                  zVectorArray *mappedNormals = nullptr,
                                  const zFloatArray *vertexHeights = nullptr,
                                  zFloatArray *mappedHeights = nullptr);

} // namespace zTsUtils

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
// clang-format off
#include <source/zToolsets/slicer/zTsUtils.cpp>
// clang-format on
#endif

#endif // ZSPACE_TS_SLICER_UTILS_H
