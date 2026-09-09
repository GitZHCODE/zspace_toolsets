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

#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <vector>

namespace zSpace {

//--------------------------
//---- PIPELINE EXECUTION METHODS
//--------------------------

// Runs m_Slicer's geometry pipeline, extracting 2D print paths.
ZSPACE_INLINE bool zTs3DP::computeSlicer() {
  m_2DPrintPaths.clear();
  m_3DPrintPaths.clear();
  bool success = m_Slicer.computeSlicer(m_2DPrintPaths);
  m_Status = m_Slicer.getStatus();
  return success;
}

// Cleans up/resamples m_2DPrintPaths in place via m_PrintSyn.
ZSPACE_INLINE bool zTs3DP::computePrintSyn() {
  bool success = m_PrintSyn.computePrintSyn(m_2DPrintPaths);
  m_Status = m_PrintSyn.getStatus();
  return success;
}

// Projects m_2DPrintPaths back onto the 3D section meshes.
ZSPACE_INLINE bool zTs3DP::computeProjection() {
  bool success = m_Slicer.computeProjection(m_2DPrintPaths, m_3DPrintPaths);
  m_Status = m_Slicer.getStatus();
  return success;
}

//--------------------------
//---- EXPORT METHOD
//--------------------------

// Exports each layer's section mesh and 3D toolpath graph to directory in OBJ format.
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
    if (zTsUtils::exportColoredMeshToOBJ(m_Slicer.getScalarField(), layerColors,
                                         path, "width_field_layer_"))
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
    if (zTsUtils::exportColoredMeshToOBJ(m_Slicer.getScalarField(), layerColors,
                                         path, "path_field_layer_"))
      exportedFiles.push_back("path_fields.obj");
  }

  // 7. 2D Toolpaths (If generated)
  if (!m_2DPrintPaths.empty()) {
    zObjGraphArray contourGraphs2D;
    contourGraphs2D.reserve(m_2DPrintPaths.size());
    for (const auto &p : m_2DPrintPaths)
      contourGraphs2D.push_back(p.graph);

    std::string path = dir + "/2D_contour_graphs.obj";
    if (zTsUtils::exportGraphToOBJ(contourGraphs2D, path, "2D_contour_layer_"))
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
