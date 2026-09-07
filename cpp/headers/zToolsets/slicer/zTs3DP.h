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
  // Pipeline Execution Methods
  //----------------------------------------------------------------------
  inline bool computeBasicStructure() {
    m_Slicer.setSettings(m_Settings);
    bool success = m_Slicer.computeBasicStructure();
    m_Status = m_Slicer.getStatus();
    return success;
  }

  inline bool computeSlice() {
    m_Slicer.setSettings(m_Settings);
    bool success = m_Slicer.computeSlice();
    m_Status = m_Slicer.getStatus();
    return success;
  }

  inline bool computeUnroll() {
    m_Slicer.setSettings(m_Settings);
    bool success = m_Slicer.computeUnroll();
    m_Status = m_Slicer.getStatus();
    return success;
  }

  inline bool computeField() {
    m_Slicer.setSettings(m_Settings);
    bool success = m_Slicer.computeField();
    m_Status = m_Slicer.getStatus();
    if (success) {
      m_Slicer.extract2DPrintPaths(m_2DPrintPaths);
    }
    return success;
  }

  inline bool computePostProcess() {
    if (m_2DPrintPaths.empty()) {
      m_Slicer.extract2DPrintPaths(m_2DPrintPaths);
    }
    // Modify PrintPaths through zTsPrintSyn without linking slicer directly
    bool success = m_PrintSyn.processAllPaths(
        m_2DPrintPaths, m_PostProcessSettings.sampleLength,
        m_PostProcessSettings.featureAngleThreshold);
    m_Status = m_PrintSyn.getStatus();
    if (!success)
      return false;

    // Slicer re-maps the modified 2D PrintPaths back to 3D section meshes
    return computePath();
  }

  inline bool computePath() {
    m_Slicer.setSettings(m_Settings);
    bool success =
        m_Slicer.project2DPrintPathsTo3D(m_2DPrintPaths, m_3DPrintPaths);
    m_Status = m_Slicer.getStatus();
    return success;
  }

  /*! \brief Runs the complete 3D printing pipeline end-to-end. */
  inline bool computePipeline() {
    if (!computeBasicStructure())
      return false;
    if (!computeSlice())
      return false;
    if (!computeUnroll())
      return false;
    if (!computeField())
      return false;
    if (!computePostProcess())
      return false;
    if (!computePath())
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

  inline zObjMesh &getUMesh() { return m_Slicer.getUMesh(); }
  inline const zObjMesh &getUMesh() const { return m_Slicer.getUMesh(); }
  inline zObjGraph &getULoop() { return m_Slicer.getULoop(); }
  inline const zObjGraph &getULoop() const { return m_Slicer.getULoop(); }

  inline zObjGraphArray &getVLoops() { return m_Slicer.getVLoops(); }
  inline const zObjGraphArray &getVLoops() const { return m_Slicer.getVLoops(); }

  inline zObjMeshArray &getSectionMeshes() {
    return m_Slicer.getSectionMeshes();
  }
  inline const zObjMeshArray &getSectionMeshes() const {
    return m_Slicer.getSectionMeshes();
  }
  inline zObjGraphArray &getSectionGraphs() {
    return m_Slicer.getSectionGraphs();
  }
  inline const zObjGraphArray &getSectionGraphs() const {
    return m_Slicer.getSectionGraphs();
  }
  inline zObjGraphArray &getBracingGraphs() {
    return m_Slicer.getBracingGraphs();
  }
  inline const zObjGraphArray &getBracingGraphs() const {
    return m_Slicer.getBracingGraphs();
  }

  inline zObjMeshArray &get2DSectionMeshes() {
    return m_Slicer.get2DSectionMeshes();
  }
  inline const zObjMeshArray &get2DSectionMeshes() const {
    return m_Slicer.get2DSectionMeshes();
  }
  inline zObjGraphArray &get2DSectionGraphs() {
    return m_Slicer.get2DSectionGraphs();
  }
  inline const zObjGraphArray &get2DSectionGraphs() const {
    return m_Slicer.get2DSectionGraphs();
  }
  inline zObjGraphArray &get2DBracingGraphs() {
    return m_Slicer.get2DBracingGraphs();
  }
  inline const zObjGraphArray &get2DBracingGraphs() const {
    return m_Slicer.get2DBracingGraphs();
  }

  inline zObjMeshScalarField &getScalarField() {
    return m_Slicer.getScalarField();
  }
  inline const zObjMeshScalarField &getScalarField() const {
    return m_Slicer.getScalarField();
  }

  inline std::vector<zScalarArray> &getPathFields() {
    return m_Slicer.getPathFields();
  }
  inline const std::vector<zScalarArray> &getPathFields() const {
    return m_Slicer.getPathFields();
  }
  inline std::vector<zScalarArray> &getWidthFields() {
    return m_Slicer.getWidthFields();
  }
  inline const std::vector<zScalarArray> &getWidthFields() const {
    return m_Slicer.getWidthFields();
  }
  inline zObjGraphArray &getTrimSlotGraphs() {
    return m_Slicer.getTrimSlotGraphs();
  }
  inline const zObjGraphArray &getTrimSlotGraphs() const {
    return m_Slicer.getTrimSlotGraphs();
  }

  // NOTE: per-layer 2D/3D toolpath data is exposed via get2DPrintPaths()/
  // get3DPrintPaths() above (this class's own m_2DPrintPaths/m_3DPrintPaths,
  // which reflect the zTsPrintSyn-cleaned result) - not by reaching into
  // getSlicer().get2DPrintPaths()/get3DPrintPaths(), which hold the slicer's
  // raw, pre-cleanup paths.

  inline const std::vector<int> &getBoundarySplitIndices() const {
    return m_Slicer.getBoundarySplitIndices();
  }
  inline const BracingGraph &getBracingParams() const {
    return m_Slicer.getBracingParams();
  }
};

/** @}*/

/** @}*/

} // namespace zSpace

#endif // ZSPACE_TS_GEOMETRY_3DP_H
