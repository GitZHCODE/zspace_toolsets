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

#ifndef ZSPACE_TS_GEOMETRY_3DP_H
#define ZSPACE_TS_GEOMETRY_3DP_H

#pragma once

#include <headers/base/zSpace_Toolsets.h>
#include <headers/zToolsets/slicer/zTsUtils.h>
#include <headers/zToolsets/slicer/zTsSDFSlicer.h>
#include <headers/zToolsets/slicer/zTsPrintSyn.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Toolsets for geometry, fabrication and computation.
 *  @{
 */

/** \addtogroup zTsGeometry
 *	\brief Geometry toolsets.
 *  @{
 */

/*! \class zTs3DP
 *	\brief Master facade and pipeline controller for robotic 3D printing toolpath generation.
 *	\details Wires together two fully decoupled engines - zTsSDFSlicer (geometry/SDF slicing)
 *           and zTsPrintSyn (toolpath cleanup) - which share no dependency on each other and
 *           only exchange data through the common PrintPath type declared in zTsUtils.h.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zTs3DP {
protected:
  //---- Input Geometry and Settings
  zObjMesh m_InputMesh;
  zObjGraph m_BracingGraph;
  SlicerSettings m_Settings;
  PostProcessSettings m_PostProcessSettings;

  //---- Sub-modules
  zTsSDFSlicer m_Slicer;
  zTsPrintSyn m_PrintSyn;

  //---- Decoupled PrintPath Containers
  std::vector<PrintPath> m_2DPrintPaths;
  std::vector<PrintPath> m_3DPrintPaths;

  //---- Status
  std::string m_Status;

public:
  //----------------------------------------------------------------------
  // Constructors & Destructors
  //----------------------------------------------------------------------
  zTs3DP() : m_Status("Initialized") {}
  ~zTs3DP() = default;

  //----------------------------------------------------------------------
  // Setters
  //----------------------------------------------------------------------
  inline void setInputMesh(zObjMesh &mesh) {
    m_InputMesh = mesh;
    m_Slicer.setInputMesh(mesh);
    m_Status = "Input mesh updated";
  }

  inline void setBracingGraph(zObjGraph &graph) {
    m_BracingGraph = graph;
    m_Slicer.setBracingGraph(graph);
    m_Status = "Bracing graph updated";
  }

  inline void setSettings(const SlicerSettings &settings) {
    m_Settings = settings;
    m_Slicer.setSettings(settings);
    m_Status = "Settings updated";
  }

  inline void setPostProcessSettings(const PostProcessSettings &settings) {
    m_PostProcessSettings = settings;
    m_Status = "Post-process settings updated";
  }

  //----------------------------------------------------------------------
  // Pipeline Execution Methods (3 steps: slicer -> print-syn -> project)
  //----------------------------------------------------------------------

  /*! \brief Runs the full zTsSDFSlicer geometry pipeline (structure, slice,
   * unroll, field) and extracts the resulting raw 2D print paths into
   * this class's own m_2DPrintPaths.
   */
  inline bool computeSlicer() {
    m_2DPrintPaths.clear();
    m_3DPrintPaths.clear();
    m_Slicer.setSettings(m_Settings);
    bool success = m_Slicer.computeBasicStructure() && m_Slicer.computeSlice() &&
                   m_Slicer.computeUnroll() &&
                   m_Slicer.computeField(m_2DPrintPaths);
    m_Status = m_Slicer.getStatus();
    return success;
  }

  /*! \brief Cleans up/resamples m_2DPrintPaths in place via zTsPrintSyn,
   * without linking the slicer and print-syn engines directly.
   */
  inline bool computePrintSyn() {
    bool success = m_PrintSyn.processAllPaths(
        m_2DPrintPaths, m_PostProcessSettings.sampleLength,
        m_PostProcessSettings.featureAngleThreshold);
    m_Status = m_PrintSyn.getStatus();
    return success;
  }

  /*! \brief Projects the (post-processed) m_2DPrintPaths back onto the 3D
   * section meshes, producing m_3DPrintPaths.
   */
  inline bool projectTo3D() {
    m_Slicer.setSettings(m_Settings);
    bool success =
        m_Slicer.projectPathsTo3D(m_2DPrintPaths, m_3DPrintPaths);
    m_Status = m_Slicer.getStatus();
    return success;
  }

  /*! \brief Runs the complete 3D printing pipeline end-to-end
   * (computeSlicer -> computePrintSyn -> projectTo3D).
   */
  inline bool computePipeline() {
    if (!computeSlicer())
      return false;
    if (!computePrintSyn())
      return false;
    if (!projectTo3D())
      return false;
    return true;
  }

  //----------------------------------------------------------------------
  // Exporters & Render Helpers
  //----------------------------------------------------------------------
  // NOTE: drawing/exporting reads print-path data that only zTs3DP owns
  // (m_2DPrintPaths/m_3DPrintPaths, the zTsPrintSyn-cleaned result), so both
  // methods live here rather than as zTsSDFSlicer members - zTsSDFSlicer no
  // longer stores print paths itself. Section meshes, V-loops, bracing
  // graphs, and the scalar field are read via m_Slicer's own getters.

  /*! \brief Builds 3D extruded rectangular pipe mesh for a specific toolpath
   * layer.
   *	\param[in]			layerIdx		- layer index.
   *	\param[out]			outMesh			- generated extruded pipe mesh.
   */
  inline void drawPathMesh(int layerIdx, zObjMesh &outMesh) const {
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
                    : m_Settings.layerSpacing;
      if (h <= 1e-4f)
        h = m_Settings.layerSpacing;

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

  /*! \brief Exports each layer's section mesh and 3D toolpath graph to
   * directory in OBJ format.
   *	\param[in]			directoryPath	- export folder path
   * (defaults to Desktop/SDF_Slicer_Export if empty).
   *	\return				bool - true if successful.
   */
  inline bool exportToOBJ(const std::string &directoryPath = "") {
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
    out << "Exported " << exportedFiles.size() << " files to " << dir
        << ": [";
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

  //----------------------------------------------------------------------
  // Getters & Data Accessors
  //----------------------------------------------------------------------
  inline bool isLoaded() const { return m_Slicer.isLoaded(); }
  inline const std::string &getStatus() const { return m_Status; }

  inline SlicerSettings &getSettings() { return m_Settings; }
  inline const SlicerSettings &getSettings() const { return m_Settings; }

  inline PostProcessSettings &getPostProcessSettings() {
    return m_PostProcessSettings;
  }
  inline const PostProcessSettings &getPostProcessSettings() const {
    return m_PostProcessSettings;
  }

  inline zTsSDFSlicer &getSlicer() { return m_Slicer; }
  inline const zTsSDFSlicer &getSlicer() const { return m_Slicer; }

  inline zTsPrintSyn &getPrintSyn() { return m_PrintSyn; }
  inline const zTsPrintSyn &getPrintSyn() const { return m_PrintSyn; }

  inline std::vector<PrintPath> &get2DPrintPaths() { return m_2DPrintPaths; }
  inline const std::vector<PrintPath> &get2DPrintPaths() const {
    return m_2DPrintPaths;
  }
  inline std::vector<PrintPath> &get3DPrintPaths() { return m_3DPrintPaths; }
  inline const std::vector<PrintPath> &get3DPrintPaths() const {
    return m_3DPrintPaths;
  }

  inline zObjMesh &getInputMesh() { return m_Slicer.getInputMesh(); }
  inline const zObjMesh &getInputMesh() const {
    return m_Slicer.getInputMesh();
  }

  inline zObjGraph &getBracingGraph() { return m_BracingGraph; }
  inline const zObjGraph &getBracingGraph() const { return m_BracingGraph; }

  // NOTE: geometry/field outputs (U-mesh, V-loops, section meshes/graphs,
  // 2D-unrolled counterparts, scalar field, path/width fields, trim slot
  // graphs, boundary split indices, bracing params) are not forwarded here -
  // access them via getSlicer().getXxx() instead, to avoid maintaining the
  // same set of getters in two places.

  // NOTE: per-layer 2D/3D toolpath data is exposed via get2DPrintPaths()/
  // get3DPrintPaths() above. zTsSDFSlicer does not store print paths itself
  // (computeField()/projectPathsTo3D() write directly into this class's own
  // m_2DPrintPaths/m_3DPrintPaths via out-parameters), so this is the only
  // place print-path data lives.

};

/** @}*/

/** @}*/

} // namespace zSpace

#endif // ZSPACE_TS_GEOMETRY_3DP_H
