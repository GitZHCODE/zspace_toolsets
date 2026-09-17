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

// Runs m_Slicer's geometry pipeline (basic structure -> slice -> unroll ->
// field).
ZSPACE_INLINE bool zTs3DP::computeSlicer() {
  m_2DPrintPaths.clear();
  m_3DPrintPaths.clear();

  if (!m_Slicer.computeBasicStructure()) {
    m_Status = m_Slicer.getStatus();
    return false;
  }
  if (!m_Slicer.computeSlice(m_SectionMeshes)) {
    m_Status = m_Slicer.getStatus();
    return false;
  }
  if (!m_Slicer.computeUnroll(m_SectionMeshes, m_SectionMeshes2D)) {
    m_Status = m_Slicer.getStatus();
    return false;
  }
  bool success = m_Slicer.computeField(m_SectionMeshes2D, m_ScalarField,
                                       m_PathFields, m_WidthFields,
                                       m_TrackFields);
  m_Status = m_Slicer.getStatus();
  return success;
}

// Runs m_PrintSyn's full toolpath pipeline (extractPrintPaths ->
// resamplePrintPath -> alignStartPoints -> projectPrintPath ->
// linkPrintPath) over the slicer's field outputs.
ZSPACE_INLINE bool zTs3DP::computePrintSyn() {
  if (!m_PrintSyn.extractPrintPaths(m_ScalarField, m_PathFields, m_WidthFields,
                                    m_TrackFields, m_2DPrintPaths)) {
    m_Status = m_PrintSyn.getStatus();
    return false;
  }
  if (!m_PrintSyn.resamplePrintPath(m_2DPrintPaths)) {
    m_Status = m_PrintSyn.getStatus();
    return false;
  }
  if (!m_PrintSyn.alignStartPoints(m_2DPrintPaths, m_SectionMeshes2D)) {
    m_Status = m_PrintSyn.getStatus();
    return false;
  }
  if (!m_PrintSyn.projectPrintPath(
          m_2DPrintPaths, m_SectionMeshes2D, m_SectionMeshes,
          m_Slicer.getVScalar(), m_Slicer.getSettings().LayerSpacing,
          m_3DPrintPaths)) {
    m_Status = m_PrintSyn.getStatus();
    return false;
  }
  bool success =
      m_PrintSyn.linkPrintPath(m_3DPrintPaths, m_ResultantPrintPath);
  m_Status = m_PrintSyn.getStatus();
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
  if (!m_SectionMeshes.empty()) {
    std::string path = dir + "/section_meshes.obj";
    if (zTsUtils::exportMeshToOBJ(m_SectionMeshes, path))
      exportedFiles.push_back("section_meshes.obj");
  }

  // 2. 3D Toolpath Graphs
  if (!m_3DPrintPaths.empty()) {
    zObjGraphArray toolpathGraphs;
    for (const auto &layerPaths : m_3DPrintPaths)
      for (const auto &p : layerPaths)
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
  if (!m_WidthFields.empty() && zFnMesh(m_ScalarField).numVertices() > 0) {
    std::vector<std::vector<zColor>> layerColors(m_WidthFields.size());
    for (size_t l = 0; l < m_WidthFields.size(); l++)
      zTsUtils::colorWidthField(m_WidthFields[l], layerColors[l]);

    std::string path = dir + "/width_fields.obj";
    if (zTsUtils::exportColoredMeshToOBJ(m_ScalarField, layerColors, path,
                                         "width_field_layer_"))
      exportedFiles.push_back("width_fields.obj");
  }

  // 6. 2D Colored Path Fields
  if (!m_PathFields.empty() && zFnMesh(m_ScalarField).numVertices() > 0) {
    std::vector<std::vector<zColor>> layerColors(m_PathFields.size());
    for (size_t l = 0; l < m_PathFields.size(); l++)
      zTsUtils::colorPathField(m_PathFields[l], layerColors[l]);

    std::string path = dir + "/path_fields.obj";
    if (zTsUtils::exportColoredMeshToOBJ(m_ScalarField, layerColors, path,
                                         "path_field_layer_"))
      exportedFiles.push_back("path_fields.obj");
  }

  // 7. 2D Toolpaths (If generated)
  if (!m_2DPrintPaths.empty()) {
    zObjGraphArray contourGraphs2D;
    for (const auto &layerPaths : m_2DPrintPaths)
      for (const auto &p : layerPaths)
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
