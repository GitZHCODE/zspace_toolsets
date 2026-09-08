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

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/slicer/zTs3DP.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <vector>

namespace zSpace {

//---- PIPELINE EXECUTION METHODS

ZSPACE_INLINE bool zTs3DP::computeSlicer() {
  m_2DPrintPaths.clear();
  m_3DPrintPaths.clear();
  bool success = m_Slicer.computeBasicStructure() && m_Slicer.computeSlice() &&
                 m_Slicer.computeUnroll() &&
                 m_Slicer.computeField(m_2DPrintPaths);
  m_Status = m_Slicer.getStatus();
  return success;
}

ZSPACE_INLINE bool zTs3DP::computePrintSyn() {
  bool success = m_PrintSyn.processAllPaths(m_2DPrintPaths);
  m_Status = m_PrintSyn.getStatus();
  return success;
}

ZSPACE_INLINE bool zTs3DP::projectTo3D() {
  bool success = m_Slicer.projectPathsTo3D(m_2DPrintPaths, m_3DPrintPaths);
  m_Status = m_Slicer.getStatus();
  return success;
}

ZSPACE_INLINE bool zTs3DP::computePipeline() {
  if (!computeSlicer())
    return false;
  if (!computePrintSyn())
    return false;
  if (!projectTo3D())
    return false;
  return true;
}

//---- DRAW / PREVIEW MESH GENERATION

ZSPACE_INLINE void zTs3DP::drawPathMesh(int layerIdx, zObjMesh &outMesh) const {
  zFnMesh fnOut(outMesh);
  fnOut.clear();

  if (layerIdx < 0 || layerIdx >= int(m_3DPrintPaths.size()))
    return;

  const PrintPath &path = m_3DPrintPaths[layerIdx];
  zFnGraph fnContour(const_cast<zObjGraph &>(path.graph));
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

  const zFloatArray *layerWeights = &path.widths;
  const zFloatArray *layerHeights = &path.heights;
  const zVectorArray *layerNormals = &path.normals;

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
                  : m_Slicer.getSettings().layerSpacing;
    if (h <= 1e-4f)
      h = m_Slicer.getSettings().layerSpacing;

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

ZSPACE_INLINE bool zTs3DP::exportToOBJ(const std::string &directoryPath) {
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
  if (!m_Slicer.getSectionMeshes().empty()) {
    std::string path = dir + "/section_meshes.obj";
    if (zTsUtils::exportMeshToOBJ(m_Slicer.getSectionMeshes(), path))
      exportedFiles.push_back("section_meshes.obj");
  }

  // 2. 3D Toolpath Graphs
  if (!m_3DPrintPaths.empty()) {
    zObjGraphArray toolpathGraphs;
    toolpathGraphs.reserve(m_3DPrintPaths.size());
    for (const auto &p : m_3DPrintPaths)
      toolpathGraphs.push_back(p.graph);

    std::string path = dir + "/toolpath_graphs.obj";
    if (zTsUtils::exportGraphToOBJ(toolpathGraphs, path, "toolpath_layer_"))
      exportedFiles.push_back("toolpath_graphs.obj");
  }

  // 3. Longitudinal V-Loops
  if (!m_Slicer.getVLoops().empty()) {
    std::string path = dir + "/vloop_graphs.obj";
    if (zTsUtils::exportGraphToOBJ(m_Slicer.getVLoops(), path, "vloop_"))
      exportedFiles.push_back("vloop_graphs.obj");
  }

  // 4. 3D Bracing Graphs
  if (!m_Slicer.getBracingGraphs().empty()) {
    std::string path = dir + "/bracing_graphs.obj";
    if (zTsUtils::exportGraphToOBJ(m_Slicer.getBracingGraphs(), path,
                                   "bracing_layer_"))
      exportedFiles.push_back("bracing_graphs.obj");
  }

  // 5. 2D Colored Width Fields (Discrete grouped colors)
  if (!m_Slicer.getWidthFields().empty() &&
      zFnMesh(m_Slicer.getScalarField()).numVertices() > 0) {
    const auto &widthFields = m_Slicer.getWidthFields();
    std::vector<std::vector<zColor>> layerColors(widthFields.size());
    for (size_t l = 0; l < widthFields.size(); l++)
      zTsUtils::colorWidthField(widthFields[l], layerColors[l]);

    std::string path = dir + "/width_fields.obj";
    if (zTsUtils::exportColoredMeshToOBJ(m_Slicer.getScalarField(),
                                         layerColors, path,
                                         "width_field_layer_"))
      exportedFiles.push_back("width_fields.obj");
  }

  // 6. 2D Colored Path Fields
  if (!m_Slicer.getPathFields().empty() &&
      zFnMesh(m_Slicer.getScalarField()).numVertices() > 0) {
    const auto &pathFields = m_Slicer.getPathFields();
    std::vector<std::vector<zColor>> layerColors(pathFields.size());
    for (size_t l = 0; l < pathFields.size(); l++)
      zTsUtils::colorPathField(pathFields[l], layerColors[l]);

    std::string path = dir + "/path_fields.obj";
    if (zTsUtils::exportColoredMeshToOBJ(m_Slicer.getScalarField(),
                                         layerColors, path,
                                         "path_field_layer_"))
      exportedFiles.push_back("path_fields.obj");
  }

  // 7. 2D Toolpaths (If generated)
  if (!m_2DPrintPaths.empty()) {
    zObjGraphArray contourGraphs2D;
    contourGraphs2D.reserve(m_2DPrintPaths.size());
    for (const auto &p : m_2DPrintPaths)
      contourGraphs2D.push_back(p.graph);

    std::string path = dir + "/2D_contour_graphs.obj";
    if (zTsUtils::exportGraphToOBJ(contourGraphs2D, path,
                                   "2D_contour_layer_"))
      exportedFiles.push_back("2D_contour_graphs.obj");
  }

  if (exportedFiles.empty()) {
    m_Status = "failed: No geometry available to export.";
    zTsUtils::printError("[zTs3DP::exportToOBJ] " + m_Status);
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
  zTsUtils::printSuccess("[zTs3DP::exportToOBJ] " + m_Status);
  return true;
}

} // namespace zSpace
