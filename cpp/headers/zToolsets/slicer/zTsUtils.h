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
// zObjMeshScalarFieldArray intentionally not typedef'd here - collides with
// zObjMeshField.h's own definition. Include zObjMeshField.h directly instead.

//--------------------------
//---- COMMON DATA STRUCTURES
//--------------------------
// SlicerSettings/PrintSynSetting live in zSlicer.h/zPrintSyn.h respectively.

/*! \struct ParameterGraph
 *	\brief Graph parameterized by 2D/3D mesh face barycentric coordinates.
 */
struct ParameterGraph {
  struct Vertex {
    int v0, v1, v2; // Triangle 3 vertex indices
    float u, v, w;  // Barycentric coordinate weights
  };

  std::vector<Vertex> vertices;
  zIntArray edgeConnects;

  bool empty() const { return vertices.empty() || edgeConnects.empty(); }
  void clear() {
    vertices.clear();
    edgeConnects.clear();
  }
};

/*! \struct PrintPath
 *	\brief Shared toolpath exchange type between zSlicer and zPrintSyn.
 */
struct PrintPath {
  zObjGraph graph;
  zFloatArray widths;
  zFloatArray heights;
  zVectorArray normals;
  // Whether this path is the double-layer boundary's outer track. Each
  // PrintPath is homogeneous (inner/outer are split into separate entries
  // upstream), so a single flag is enough. Single-perimeter/bracing-only
  // paths default false (treated as "inner").
  bool isOuterTrack = false;
};

//--------------------------
//---- COMMON GEOMETRY / TOPOLOGY UTILITIES
//--------------------------

namespace zTsUtils {

//---- Thread-Safe Console Logging

void printSuccess(const std::string &msg);
void printError(const std::string &msg);

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

//---- Single-Chain Graph Closed/Open Detection

/*! \brief Infers closed/open from a wraparound edge (vertex 0 <-> N-1). */
bool isClosedGraph(const zObjGraph &graph);

//---- Universal Polyline Sampling by Arc-Length with Optional Tangent Vector

bool samplePolyline(const zPoint *positions, const zScalar *scalars,
                    size_t count, double t, zPoint &outPoint,
                    zVector *outTangent = nullptr);
bool samplePolyline(const zObjGraph &polylineGraph, const zScalarArray &scalars,
                    double t, zPoint &outPoint, zVector *outTangent = nullptr);

//---- Universal Mesh Barycentric Parameterizer (2D flat and 3D space meshes)

bool parameterizeGraph(const zObjGraph &inGraph, zObjMesh &inMesh,
                             ParameterGraph &outParams);

//---- Universal Barycentric Mesh Graph Reconstructor

bool realizeParameterGraph(
    const ParameterGraph &inParams, zObjMesh &targetMesh, zObjGraph &outGraph,
    zVectorArray *outNormals = nullptr,
    const zFloatArray *vertexScalars = nullptr,
    zFloatArray *outHeights = nullptr);

//---- Barycentric Mesh-to-Mesh Graph Projection

bool remapGraph(zObjGraph &graph, zObjMesh &inMesh2D,
                                  zObjMesh &projectionMesh3D,
                                  zVectorArray *mappedNormals = nullptr,
                                  const zFloatArray *vertexHeights = nullptr,
                                  zFloatArray *mappedHeights = nullptr);

//---- Toolpath Pipe Mesh Construction

/*! \brief Builds a 3D extruded pipe mesh from a print path. */
void buildPathMesh(const PrintPath &path, float defaultLayerHeight,
                   zObjMesh &outMesh);

//--------------------------
//---- FIELD-MESH COLORING & OBJ EXPORT UTILITIES
//--------------------------

//---- Color Format Conversion

zColor hslToRgb(float h, float s, float l);

//---- Distinct Shuffled Color Palette Generator

void generateColorPalette(int numColors, std::vector<zColor> &outColors,
                          unsigned int seed = 1337, float saturation = 0.85f,
                          float lightness = 0.50f);

//---- Field-Mesh Coloring: Positive Values Solid, Negative Values a
//---- Gradient Between Two Input Colors (near-zero -> most-negative)

void colorPathField(
    const zScalarArray &values, std::vector<zColor> &outColors,
    const zColor &positiveColor = zColor(1.0f, 1.0f, 1.0f, 1.0f),
    const zColor &negativeColorNearZero = zColor(0.0f, 40.0f / 255.0f,
                                                 240.0f / 255.0f, 1.0f),
    const zColor &negativeColorAtMin = zColor(180.0f / 255.0f, 200.0f / 255.0f,
                                              1.0f, 1.0f));
void colorPathField(
    zObjMesh &mesh, const zScalarArray &values,
    const zColor &positiveColor = zColor(1.0f, 1.0f, 1.0f, 1.0f),
    const zColor &negativeColorNearZero = zColor(0.0f, 40.0f / 255.0f,
                                                 240.0f / 255.0f, 1.0f),
    const zColor &negativeColorAtMin = zColor(180.0f / 255.0f, 200.0f / 255.0f,
                                              1.0f, 1.0f));

//---- Field-Mesh Coloring: Discrete Lookup by Index Into a Generated Color
//---- List (one distinct color per unique value)

void colorWidthField(const zScalarArray &values,
                     std::vector<zColor> &outColors);
void colorWidthField(zObjMesh &mesh, const zScalarArray &values);

//---- Generic Multi-Object Wavefront OBJ Exporters

bool exportMeshToOBJ(const zObjMeshArray &meshes, const std::string &filePath);
bool exportGraphToOBJ(const zObjGraphArray &graphs, const std::string &filePath,
                      const std::string &objectPrefix = "layer_");
bool exportColoredMeshToOBJ(const zObjMesh &baseMesh,
                            const std::vector<std::vector<zColor>> &layerColors,
                            const std::string &filePath,
                            const std::string &objectPrefix = "layer_");

//---- Wavefront OBJ Graph Importer

/*! \brief Reads a graph's v/l lines from an OBJ file (all objects merged). */
bool readGraphFromOBJ(const std::string &filePath, zObjGraph &outGraph);

//---- Generic Mesh/Graph Transform & Segmentation

/*! \brief Translates all vertex positions of a mesh. */
void translateMesh(zObjMesh &mesh, const zVector &offset);

/*! \brief Translates all vertex positions of a graph. */
void translateGraph(zObjGraph &graph, const zVector &offset);

/*! \brief Segments a graph at specified split vertex indices. */
bool segmentGraph(const zObjGraph &sourceBoundaryGraph,
                  const std::vector<int> &splitIndices,
                  std::vector<zObjGraph> &outBoundarySegments);

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
#include <source/zToolsets/slicer/zTsUtils_Export.cpp>
// clang-format on
#endif

#endif // ZSPACE_TS_SLICER_UTILS_H
