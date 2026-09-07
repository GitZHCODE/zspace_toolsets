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

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <vector>

namespace zSpace {

using namespace zTsUtils;

//---- SCALAR FIELD COLORING METHODS (zTsSDFSlicer)
//---- Thin forwarders onto the generic implementations in zTsUtils (see
//---- zTsUtils_Export.cpp), kept here so existing callers of
//---- zTsSDFSlicer::colorPathField/colorWidthField/... keep working.
//---- colorPathField forwards to zTsUtils::colorField with a fixed
//---- purple -> yellow gradient (its own start/end colors are not exposed
//---- as parameters, to preserve the existing 2-argument public API).

ZSPACE_INLINE void
zTsSDFSlicer::colorPathField(const zScalarArray &values,
                             std::vector<zColor> &outColors) {
  static const zColor kStartColor(0.45f, 0.10f, 0.85f, 1.0f); // purple
  static const zColor kEndColor(1.00f, 0.85f, 0.15f, 1.0f);   // yellow
  zTsUtils::colorField(values, kStartColor, kEndColor, outColors);
}

ZSPACE_INLINE void zTsSDFSlicer::colorPathField(zObjMesh &mesh,
                                                const zScalarArray &values) {
  static const zColor kStartColor(0.45f, 0.10f, 0.85f, 1.0f); // purple
  static const zColor kEndColor(1.00f, 0.85f, 0.15f, 1.0f);   // yellow
  zTsUtils::colorField(mesh, values, kStartColor, kEndColor);
}

ZSPACE_INLINE zColor zTsSDFSlicer::hslToRgb(float h, float s, float l) {
  return zTsUtils::hslToRgb(h, s, l);
}

ZSPACE_INLINE void
zTsSDFSlicer::generateColorPalette(int numColors, std::vector<zColor> &outColors,
                                   unsigned int seed, float saturation,
                                   float lightness) {
  zTsUtils::generateColorPalette(numColors, outColors, seed, saturation,
                                 lightness);
}

ZSPACE_INLINE void
zTsSDFSlicer::colorWidthField(const zScalarArray &values,
                              std::vector<zColor> &outColors) {
  zTsUtils::colorWidthField(values, outColors);
}

ZSPACE_INLINE void zTsSDFSlicer::colorWidthField(zObjMesh &mesh,
                                                 const zScalarArray &values) {
  zTsUtils::colorWidthField(mesh, values);
}

//---- GEOMETRY HELPERS

ZSPACE_INLINE void zTsSDFSlicer::translateMesh(zObjMesh &mesh,
                                               const zVector &offset) {
  zFnMesh fn(mesh);
  if (fn.numVertices() <= 0)
    return;
  zPointArray pts;
  fn.getVertexPositions(pts);
  for (auto &p : pts) {
    p.x += offset.x;
    p.y += offset.y;
    p.z += offset.z;
  }
  fn.setVertexPositions(pts);
}

ZSPACE_INLINE void zTsSDFSlicer::translateGraph(zObjGraph &graph,
                                                const zVector &offset) {
  zFnGraph fn(graph);
  if (fn.numVertices() <= 0)
    return;
  zPointArray pts;
  fn.getVertexPositions(pts);
  for (auto &p : pts) {
    p.x += offset.x;
    p.y += offset.y;
    p.z += offset.z;
  }
  fn.setVertexPositions(pts);
}

ZSPACE_INLINE bool
zTsSDFSlicer::segmentGraph(const zObjGraph &sourceBoundaryGraph,
                           const std::vector<int> &splitIndices,
                           std::vector<zObjGraph> &outBoundarySegments) {
  outBoundarySegments.clear();
  zFnGraph fnSource(const_cast<zObjGraph &>(sourceBoundaryGraph));
  zPointArray positions;
  fnSource.getVertexPositions(positions);

  const int N = int(positions.size());
  if (N < 2)
    return false;

  const int numSegs = int(splitIndices.size());
  if (numSegs == 0) {
    outBoundarySegments.push_back(sourceBoundaryGraph);
    return true;
  }

  outBoundarySegments.reserve(numSegs);

  for (int seg = 0; seg < numSegs; seg++) {
    const int startIdx = splitIndices[seg];
    const int endIdx =
        (seg + 1 < numSegs) ? splitIndices[seg + 1] : splitIndices[0];

    zPointArray outPositions;
    zIntArray outEdges;

    if (startIdx <= endIdx) {
      for (int k = startIdx; k <= endIdx; k++) {
        outPositions.push_back(positions[k]);
      }
    } else {
      for (int k = startIdx; k < N; k++) {
        outPositions.push_back(positions[k]);
      }
      for (int k = 0; k <= endIdx; k++) {
        outPositions.push_back(positions[k]);
      }
    }

    if (outPositions.size() >= 2) {
      for (size_t e = 0; e + 1 < outPositions.size(); e++) {
        outEdges.push_back(int(e));
        outEdges.push_back(int(e + 1));
      }

      zObjGraph segGraph;
      zFnGraph fnSeg(segGraph);
      fnSeg.create(outPositions, outEdges);
      outBoundarySegments.push_back(segGraph);
    }
  }

  return !outBoundarySegments.empty();
}

ZSPACE_INLINE void
zTsSDFSlicer::drawPathMesh(int layerIdx, zObjMesh &outMesh) const {
  zFnMesh fnOut(outMesh);
  fnOut.clear();

  if (layerIdx < 0 || layerIdx >= int(contourGraphs.size()))
    return;

  zFnGraph fnContour(const_cast<zObjGraph &>(contourGraphs[layerIdx]));
  zPointArray points;
  zIntArray edges;
  fnContour.getVertexPositions(points);
  fnContour.getEdgeData(edges);
  const int N = int(points.size());
  if (N < 2)
    return;

  bool closed = false;
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    if ((edges[e] == 0 && edges[e + 1] == N - 1) ||
        (edges[e + 1] == 0 && edges[e] == N - 1)) {
      closed = true;
      break;
    }
  }

  const zFloatArray *layerWeights =
      (layerIdx < int(contourWeights.size())) ? &contourWeights[layerIdx]
                                             : nullptr;
  const zFloatArray *layerHeights =
      (layerIdx < int(contourHeights.size())) ? &contourHeights[layerIdx]
                                             : nullptr;
  const zVectorArray *layerNormals =
      (layerIdx < int(contourNormals.size())) ? &contourNormals[layerIdx]
                                             : nullptr;

  zPointArray meshPositions;
  zIntArray polyCounts;
  zIntArray polyConnects;
  meshPositions.reserve(N * 4);

  for (int i = 0; i < N; i++) {
    zPoint p = points[i];

    // 1. Tangent direction
    zVector pathDir;
    if (closed) {
      int prevIdx = (i - 1 + N) % N;
      int nextIdx = (i + 1) % N;
      pathDir = zVector(points[nextIdx].x - points[prevIdx].x,
                        points[nextIdx].y - points[prevIdx].y,
                        points[nextIdx].z - points[prevIdx].z);
    } else {
      if (i == 0) {
        pathDir = zVector(points[1].x - points[0].x,
                          points[1].y - points[0].y,
                          points[1].z - points[0].z);
      } else if (i == N - 1) {
        pathDir = zVector(points[N - 1].x - points[N - 2].x,
                          points[N - 1].y - points[N - 2].y,
                          points[N - 1].z - points[N - 2].z);
      } else {
        pathDir = zVector(points[i + 1].x - points[i - 1].x,
                          points[i + 1].y - points[i - 1].y,
                          points[i + 1].z - points[i - 1].z);
      }
    }
    if (pathDir.length() < 1e-6)
      pathDir = zVector(1, 0, 0);
    pathDir.normalize();

    // 2. 3D surface normal
    zVector norm = (layerNormals && i < int(layerNormals->size()))
                       ? (*layerNormals)[i]
                       : zVector(0, 0, 1);
    if (norm.length() < 1e-6)
      norm = zVector(0, 0, 1);
    norm.normalize();

    // 3. In-plane binormal perpendicular to tangent
    zVector binormal = norm ^ pathDir;
    if (binormal.length() < 1e-6) {
      zVector temp(0, 0, 1);
      if (std::abs(pathDir.z) > 0.9)
        temp = zVector(0, 1, 0);
      binormal = temp ^ pathDir;
    }
    binormal.normalize();

    // 4. Per-vertex width and height
    float w = (layerWeights && i < int(layerWeights->size()))
                  ? (*layerWeights)[i]
                  : 0.20f;
    if (w <= 1e-4f)
      w = 0.20f;
    float h = (layerHeights && i < int(layerHeights->size()))
                  ? (*layerHeights)[i]
                  : i_Settings.layerSpacing;
    if (h <= 1e-4f)
      h = i_Settings.layerSpacing;

    // 5. Generate 4 corner vertices for rectangular cross section
    const float halfW = w * 0.5f;
    zPoint c0(p.x + binormal.x * halfW, p.y + binormal.y * halfW,
              p.z + binormal.z * halfW); // Top-Right
    zPoint c1(p.x - binormal.x * halfW, p.y - binormal.y * halfW,
              p.z - binormal.z * halfW); // Top-Left
    zPoint c2(c1.x - norm.x * h, c1.y - norm.y * h,
              c1.z - norm.z * h); // Bottom-Left
    zPoint c3(c0.x - norm.x * h, c0.y - norm.y * h,
              c0.z - norm.z * h); // Bottom-Right

    meshPositions.push_back(c0);
    meshPositions.push_back(c1);
    meshPositions.push_back(c2);
    meshPositions.push_back(c3);
  }

  // 6. Construct 4 quads per segment along the toolpath tube
  const int segCount = closed ? N : (N - 1);
  for (int i = 0; i < segCount; i++) {
    int next = (i + 1) % N;
    // Top face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 0);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(i * 4 + 1);

    // Left face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 1);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(i * 4 + 2);

    // Bottom face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 2);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(i * 4 + 3);

    // Right face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 3);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(i * 4 + 0);
  }

  // Cap ends for open curves
  if (!closed && N >= 2) {
    polyCounts.push_back(4);
    polyConnects.push_back(0);
    polyConnects.push_back(3);
    polyConnects.push_back(2);
    polyConnects.push_back(1);

    int last = (N - 1) * 4;
    polyCounts.push_back(4);
    polyConnects.push_back(last + 0);
    polyConnects.push_back(last + 1);
    polyConnects.push_back(last + 2);
    polyConnects.push_back(last + 3);
  }

  fnOut.create(meshPositions, polyCounts, polyConnects);
}

//---- EXPORT METHOD

ZSPACE_INLINE bool zTsSDFSlicer::exportToOBJ(const std::string &directoryPath) {
  std::string dir = directoryPath;
  if (dir.empty()) {
    const char *userProfile = std::getenv("USERPROFILE");
    dir = userProfile
              ? (std::string(userProfile) + "/Desktop/SDF_Slicer_Export")
              : "C:/Users/Default/Desktop/SDF_Slicer_Export";
  }

  std::error_code ec;
  std::filesystem::create_directories(dir, ec);

  std::vector<std::string> exportedFiles;

  // 1. 3D Section Meshes
  if (!sectionMeshes.empty()) {
    std::string path = dir + "/section_meshes.obj";
    if (zTsUtils::exportMeshesToMultiOBJ(sectionMeshes, path))
      exportedFiles.push_back("section_meshes.obj");
  }

  // 2. 3D Toolpath Graphs
  if (!contourGraphs.empty()) {
    std::string path = dir + "/toolpath_graphs.obj";
    if (zTsUtils::exportGraphsToMultiOBJ(contourGraphs, path, "toolpath_layer_"))
      exportedFiles.push_back("toolpath_graphs.obj");
  }

  // 3. Longitudinal V-Loops
  if (!m_VLoops.empty()) {
    std::string path = dir + "/vloop_graphs.obj";
    if (zTsUtils::exportGraphsToMultiOBJ(m_VLoops, path, "vloop_"))
      exportedFiles.push_back("vloop_graphs.obj");
  }

  // 4. 3D Bracing Graphs
  if (!bracingGraphs.empty()) {
    std::string path = dir + "/bracing_graphs.obj";
    if (zTsUtils::exportGraphsToMultiOBJ(bracingGraphs, path, "bracing_layer_"))
      exportedFiles.push_back("bracing_graphs.obj");
  }

  // 5. 2D Colored Width Fields (Discrete grouped colors)
  if (!widthFields.empty() && zFnMesh(m_ScalarField).numVertices() > 0) {
    std::vector<std::vector<zColor>> layerColors(widthFields.size());
    for (size_t l = 0; l < widthFields.size(); l++)
      zTsUtils::colorWidthField(widthFields[l], layerColors[l]);

    std::string path = dir + "/width_fields.obj";
    if (zTsUtils::exportColoredMeshesToMultiOBJ(m_ScalarField, layerColors, path,
                                                "width_field_layer_"))
      exportedFiles.push_back("width_fields.obj");
  }

  // 6. 2D Colored Path Fields
  if (!pathFields.empty() && zFnMesh(m_ScalarField).numVertices() > 0) {
    const zColor pathFieldStartColor(0.45f, 0.10f, 0.85f, 1.0f); // purple
    const zColor pathFieldEndColor(1.00f, 0.85f, 0.15f, 1.0f);   // yellow
    std::vector<std::vector<zColor>> layerColors(pathFields.size());
    for (size_t l = 0; l < pathFields.size(); l++)
      zTsUtils::colorField(pathFields[l], pathFieldStartColor,
                           pathFieldEndColor, layerColors[l]);

    std::string path = dir + "/path_fields.obj";
    if (zTsUtils::exportColoredMeshesToMultiOBJ(m_ScalarField, layerColors, path,
                                                "path_field_layer_"))
      exportedFiles.push_back("path_fields.obj");
  }

  // 7. 2D Toolpaths (If generated)
  if (!_2DContourGraphs.empty()) {
    std::string path = dir + "/2D_contour_graphs.obj";
    if (zTsUtils::exportGraphsToMultiOBJ(_2DContourGraphs, path,
                                         "2D_contour_layer_"))
      exportedFiles.push_back("2D_contour_graphs.obj");
  }

  if (exportedFiles.empty()) {
    m_Status = "failed: No geometry available to export.";
    printError("[exportToOBJ] " + m_Status);
    return false;
  }

  std::ostringstream out;
  out << "Exported " << exportedFiles.size() << " files to " << dir << ": [";
  for (size_t i = 0; i < exportedFiles.size(); i++) {
    if (i > 0)
      out << ", ";
    out << exportedFiles[i];
  }
  out << "]";
  m_Status = out.str();
  printSuccess("[exportToOBJ] " + m_Status);
  return true;
}

} // namespace zSpace
