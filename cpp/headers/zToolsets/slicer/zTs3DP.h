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
#include <headers/zToolsets/slicer/zPrintSyn.h>
#include <headers/zToolsets/slicer/zSlicer.h>
#include <headers/zToolsets/slicer/zTsUtils.h>

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
 *	\brief Facade wiring zSlicer (SDF slicing) and zPrintSyn (toolpath
 * cleanup) into one pipeline.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zTs3DP {
protected:
  zObjMesh m_InputMesh;
  zObjGraph m_BracingGraph;
  zSlicer m_Slicer;
  zPrintSyn m_PrintSyn;
  std::vector<std::vector<PrintPath>> m_2DPrintPaths;
  std::vector<std::vector<PrintPath>> m_3DPrintPaths;
  PrintPath m_ResultantPrintPath;

  zObjMeshArray m_SectionMeshes;
  zObjMeshArray m_SectionMeshes2D;
  zObjMeshScalarField m_ScalarField;
  std::vector<zScalarArray> m_PathFields;
  std::vector<zScalarArray> m_WidthFields;
  std::vector<zScalarArray> m_TrackFields;

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

  // getSettings()/getPrintSynSetting() below return mutable references into
  // m_Slicer/m_PrintSyn - assign through those instead of a setter.

  //----------------------------------------------------------------------
  // Pipeline Execution Methods
  //----------------------------------------------------------------------

  /*! \brief Runs m_Slicer's geometry pipeline (basic structure -> slice ->
   * unroll -> field), producing m_SectionMeshes/m_SectionMeshes2D/
   * m_ScalarField/m_PathFields/m_WidthFields/m_TrackFields. */
  bool computeSlicer();

  /*! \brief Runs m_PrintSyn's full toolpath pipeline (extractPrintPaths ->
   * resamplePrintPath -> alignStartPoints -> projectPrintPath ->
   * linkPrintPath) over the slicer's field outputs, producing
   * m_2DPrintPaths/m_3DPrintPaths/m_ResultantPrintPath. Call after
   * computeSlicer. */
  bool computePrintSyn();

  /*! \brief Runs computeSlicer -> computePrintSyn. */
  bool computePipeline() {
    if (!computeSlicer())
      return false;
    if (!computePrintSyn())
      return false;
    return true;
  }

  //----------------------------------------------------------------------
  // Exporters
  //----------------------------------------------------------------------

  /*! \brief Exports each layer's section mesh and 3D toolpath graph to
   * directory in OBJ format. */
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

  inline zSlicer &getSlicer() { return m_Slicer; }
  inline const zSlicer &getSlicer() const { return m_Slicer; }

  inline zPrintSyn &getPrintSyn() { return m_PrintSyn; }
  inline const zPrintSyn &getPrintSyn() const { return m_PrintSyn; }

  inline std::vector<std::vector<PrintPath>> &get2DPrintPaths() {
    return m_2DPrintPaths;
  }
  inline const std::vector<std::vector<PrintPath>> &get2DPrintPaths() const {
    return m_2DPrintPaths;
  }
  inline std::vector<std::vector<PrintPath>> &get3DPrintPaths() {
    return m_3DPrintPaths;
  }
  inline const std::vector<std::vector<PrintPath>> &get3DPrintPaths() const {
    return m_3DPrintPaths;
  }
  inline PrintPath &getResultantPrintPath() { return m_ResultantPrintPath; }
  inline const PrintPath &getResultantPrintPath() const {
    return m_ResultantPrintPath;
  }

  inline zObjMeshArray &getSectionMeshes() { return m_SectionMeshes; }
  inline const zObjMeshArray &getSectionMeshes() const {
    return m_SectionMeshes;
  }
  inline zObjMeshArray &getSectionMeshes2D() { return m_SectionMeshes2D; }
  inline const zObjMeshArray &getSectionMeshes2D() const {
    return m_SectionMeshes2D;
  }
  inline zObjMeshScalarField &getScalarField() { return m_ScalarField; }
  inline const zObjMeshScalarField &getScalarField() const {
    return m_ScalarField;
  }
  inline std::vector<zScalarArray> &getPathFields() { return m_PathFields; }
  inline const std::vector<zScalarArray> &getPathFields() const {
    return m_PathFields;
  }
  inline std::vector<zScalarArray> &getWidthFields() { return m_WidthFields; }
  inline const std::vector<zScalarArray> &getWidthFields() const {
    return m_WidthFields;
  }
  inline std::vector<zScalarArray> &getTrackFields() { return m_TrackFields; }
  inline const std::vector<zScalarArray> &getTrackFields() const {
    return m_TrackFields;
  }

  inline zObjMesh &getInputMesh() { return m_Slicer.getInputMesh(); }
  inline const zObjMesh &getInputMesh() const {
    return m_Slicer.getInputMesh();
  }

  inline zObjGraph &getBracingGraph() { return m_BracingGraph; }
  inline const zObjGraph &getBracingGraph() const { return m_BracingGraph; }

  // Other geometry/field outputs (U-mesh, V-loops, section meshes/graphs,
  // scalar field, etc.) are accessed via getSlicer().getXxx() instead of
  // being forwarded here.
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
