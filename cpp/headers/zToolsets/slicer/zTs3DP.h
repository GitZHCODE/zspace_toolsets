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
    m_Slicer.setSettings(settings);
    m_Status = "Settings updated";
  }

  inline void setPrintSynSetting(const PrintSynSetting &settings) {
    m_PrintSyn.setSettings(settings);
    m_Status = "Post-process settings updated";
  }

  //----------------------------------------------------------------------
  // Pipeline Execution Methods (3 steps: slicer -> print-syn -> project)
  //----------------------------------------------------------------------

  /*! \brief Runs the full zTsSDFSlicer geometry pipeline (structure, slice,
   * unroll, field) and extracts the resulting raw 2D print paths into
   * this class's own m_2DPrintPaths.
   */
  bool computeSlicer();

  /*! \brief Cleans up/resamples m_2DPrintPaths in place via zTsPrintSyn,
   * without linking the slicer and print-syn engines directly.
   */
  bool computePrintSyn();

  /*! \brief Projects the (post-processed) m_2DPrintPaths back onto the 3D
   * section meshes, producing m_3DPrintPaths.
   */
  bool projectTo3D();

  /*! \brief Runs the complete 3D printing pipeline end-to-end
   * (computeSlicer -> computePrintSyn -> projectTo3D).
   */
  bool computePipeline();

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
  void drawPathMesh(int layerIdx, zObjMesh &outMesh) const;

  /*! \brief Exports each layer's section mesh and 3D toolpath graph to
   * directory in OBJ format.
   *	\param[in]			directoryPath	- export folder path
   * (defaults to Desktop/SDF_Slicer_Export if empty).
   *	\return				bool - true if successful.
   */
  bool exportToOBJ(const std::string &directoryPath = "");

  //----------------------------------------------------------------------
  // Getters & Data Accessors
  //----------------------------------------------------------------------
  inline bool isLoaded() const { return m_Slicer.isLoaded(); }
  inline const std::string &getStatus() const { return m_Status; }

  inline SlicerSettings &getSettings() { return m_Slicer.getSettings(); }
  inline const SlicerSettings &getSettings() const {
    return m_Slicer.getSettings();
  }

  inline PrintSynSetting &getPrintSynSetting() {
    return m_PrintSyn.getSettings();
  }
  inline const PrintSynSetting &getPrintSynSetting() const {
    return m_PrintSyn.getSettings();
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

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include <source/zToolsets/slicer/zTs3DP.cpp>
#endif

#endif // ZSPACE_TS_GEOMETRY_3DP_H
