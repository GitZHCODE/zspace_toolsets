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

#include <iostream>
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
    m_Slicer.setSettings(m_Settings);
    bool success = m_Slicer.computeBasicStructure() && m_Slicer.computeSlice() &&
                   m_Slicer.computeUnroll() && m_Slicer.computeField();
    m_Status = m_Slicer.getStatus();
    if (success) {
      m_2DPrintPaths = m_Slicer.get2DPrintPaths();
    }
    return success;
  }

  /*! \brief Cleans up/resamples m_2DPrintPaths in place via zTsPrintSyn,
   * without linking the slicer and print-syn engines directly.
   */
  inline bool computePrintSyn() {
    if (m_2DPrintPaths.empty()) {
      m_2DPrintPaths = m_Slicer.get2DPrintPaths();
    }
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
  inline bool exportToOBJ(const std::string &directoryPath = "") {
    return m_Slicer.exportToOBJ(directoryPath);
  }

  inline void drawPathMesh(int layerIdx, zObjMesh &outMesh) const {
    m_Slicer.drawPathMesh(layerIdx, outMesh);
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
  // get3DPrintPaths() above (this class's own m_2DPrintPaths/m_3DPrintPaths,
  // which reflect the zTsPrintSyn-cleaned result) - not by reaching into
  // getSlicer().get2DPrintPaths()/get3DPrintPaths(), which hold the slicer's
  // raw, pre-cleanup paths.

};

/** @}*/

/** @}*/

} // namespace zSpace

#endif // ZSPACE_TS_GEOMETRY_3DP_H
