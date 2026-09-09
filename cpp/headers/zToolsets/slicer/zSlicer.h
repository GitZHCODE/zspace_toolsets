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

#ifndef ZSPACE_TS_GEOMETRY_SLICER_H
#define ZSPACE_TS_GEOMETRY_SLICER_H

#pragma once

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/base/zSpace_Toolsets.h>
#include <headers/zToolsets/slicer/zTsUtils.h>

#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnMeshField.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>
#include <zspace/zInterface/objects/zObjMeshField.h>

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

/*! \struct SlicerSettings
 *	\brief Configuration settings for geodesic slicing and SDF/toolpath generation.
 */
struct SlicerSettings {
  int cornerStartVertexId = 21;
  float layerSpacing = 1.0f;
  bool topoSplit = true;
  bool angleSplit = true;
  float geometricAngleThreshold = 90.0f; // degrees (corner angle < threshold)

  // SDF & Toolpath Parameters
  bool doublelayerboundary = false;
  std::vector<float> printBoundaryWidth = {0.5f};
  float printBracingWidth = 0.5f;
  float printOverlapWidth = 0.0f;
  float trimSlotStaggerEven = 0.5f;
  float trimSlotStaggerOdd = 0.5f;
  int boundaryTrimSegmentId = 1;
  float boundaryTrimSlotRatioEven = 0.45f;
  float boundaryTrimSlotRatioOdd = 0.55f;

  float fieldBoundingScale = 1.05f;
  float fieldCellSize = 0.05f; // Target square grid cell size (default 0.05)
};

/*! \class zSlicer
 *	\brief SDF based geodesic slicing and bracing generation engine.
 */
class ZSPACE_TOOLSETS zSlicer {
protected:
  //--------------------------
  //---- PROTECTED ATTRIBUTES
  //--------------------------

  //---- INPUTS (i_)
  zObjMesh *i_InputMesh;
  zObjGraph *i_BracingGraph;
  SlicerSettings i_Settings;

  //---- MIDDLE / REUSABLE DATA (m_)
  zObjMesh m_UMesh;
  zObjGraph m_ULoop;
  zObjGraphArray m_VLoops;
  std::vector<zScalarArray> m_VScalar;
  zObjMeshScalarField m_ScalarField;
  std::vector<int> m_BoundarySplitIndices;
  ParameterGraph m_BracingParams;
  std::string m_Status;

  //---- OUTPUTS
  zObjMeshArray sectionMeshes;
  zObjGraphArray sectionGraphs;
  zObjGraphArray bracingGraphs;
  zObjMeshArray sectionMeshes2D;
  zObjGraphArray sectionGraphs2D;
  zObjGraphArray bracingGraphs2D;
  std::vector<zScalarArray> pathFields;
  std::vector<zScalarArray> widthFields;
  zObjGraphArray trimSlotGraphs;

public:
  //--------------------------
  //---- CONSTRUCTOR & DESTRUCTOR
  //--------------------------

  zSlicer() {
    i_InputMesh = nullptr;
    i_BracingGraph = nullptr;
    m_Status = "Ready";
  }

  explicit zSlicer(const SlicerSettings &inSettings) : zSlicer() {
    i_Settings = inSettings;
  }

  ~zSlicer() {}

  //--------------------------
  //---- SET METHODS
  //--------------------------

  void setSettings(const SlicerSettings &inSettings) {
    i_Settings = inSettings;
  }
  const SlicerSettings &getSettings() const { return i_Settings; }
  SlicerSettings &getSettings() { return i_Settings; }

  void setInputMesh(zObjMesh &inMesh) { i_InputMesh = &inMesh; }
  void setBracingGraph(zObjGraph &inGraph) { i_BracingGraph = &inGraph; }

  void setStatus(const std::string &inStatus) { m_Status = inStatus; }
  const std::string &getStatus() const { return m_Status; }

  //--------------------------
  //---- GET METHODS
  //--------------------------

  bool isLoaded() const { return i_InputMesh != nullptr; }
  int getLayerCount() const { return int(sectionMeshes.size()); }

  zObjMesh &getInputMesh() { return *i_InputMesh; }
  const zObjMesh &getInputMesh() const { return *i_InputMesh; }

  zObjMesh &getUMesh() { return m_UMesh; }
  const zObjMesh &getUMesh() const { return m_UMesh; }
  zObjGraph &getULoop() { return m_ULoop; }
  const zObjGraph &getULoop() const { return m_ULoop; }

  zObjGraphArray &getVLoops() { return m_VLoops; }
  const zObjGraphArray &getVLoops() const { return m_VLoops; }

  std::vector<zScalarArray> &getVScalar() { return m_VScalar; }
  const std::vector<zScalarArray> &getVScalar() const { return m_VScalar; }

  zObjMeshArray &getSectionMeshes() { return sectionMeshes; }
  const zObjMeshArray &getSectionMeshes() const { return sectionMeshes; }
  zObjGraphArray &getSectionGraphs() { return sectionGraphs; }
  const zObjGraphArray &getSectionGraphs() const { return sectionGraphs; }
  zObjGraphArray &getBracingGraphs() { return bracingGraphs; }
  const zObjGraphArray &getBracingGraphs() const { return bracingGraphs; }

  zObjMeshArray &get2DSectionMeshes() { return sectionMeshes2D; }
  const zObjMeshArray &get2DSectionMeshes() const { return sectionMeshes2D; }
  zObjGraphArray &get2DSectionGraphs() { return sectionGraphs2D; }
  const zObjGraphArray &get2DSectionGraphs() const { return sectionGraphs2D; }
  zObjGraphArray &get2DBracingGraphs() { return bracingGraphs2D; }
  const zObjGraphArray &get2DBracingGraphs() const { return bracingGraphs2D; }

  std::vector<zScalarArray> &getPathFields() { return pathFields; }
  const std::vector<zScalarArray> &getPathFields() const { return pathFields; }
  std::vector<zScalarArray> &getWidthFields() { return widthFields; }
  const std::vector<zScalarArray> &getWidthFields() const {
    return widthFields;
  }

  zObjMeshScalarField &getScalarField() { return m_ScalarField; }
  const zObjMeshScalarField &getScalarField() const { return m_ScalarField; }

  zObjGraphArray &getTrimSlotGraphs() { return trimSlotGraphs; }
  const zObjGraphArray &getTrimSlotGraphs() const { return trimSlotGraphs; }

  const std::vector<int> &getBoundarySplitIndices() const {
    return m_BoundarySplitIndices;
  }
  const ParameterGraph &getBracingParams() const { return m_BracingParams; }

  //--------------------------
  //---- MAIN COMPUTE METHODS
  //--------------------------

  /*! \brief Corner extraction, V-loops, U-mesh, U-loop, bracing parameterization. */
  bool computeBasicStructure();

  /*! \brief Calculates layer meshes, boundary section graphs, layer bracing graphs. */
  bool computeSlice();

  /*! \brief 2D unrolling for all 3D section meshes, boundary graphs, bracing graphs. */
  bool computeUnroll();

  /*! \brief 2D SDF scalar fields; extracts flat isocontour print paths into outPrintPaths2D. */
  bool computeField(std::vector<PrintPath> &outPrintPaths2D);

  /*! \brief Runs computeBasicStructure -> computeSlice -> computeUnroll -> computeField. */
  bool computeSlicer(std::vector<PrintPath> &outPrintPaths2D) {
    return computeBasicStructure() && computeSlice() && computeUnroll() &&
           computeField(outPrintPaths2D);
  }

  /*! \brief Projects 2D PrintPaths back to 3D section meshes. */
  bool projectPathsTo3D(const std::vector<PrintPath> &in2DPaths,
                        std::vector<PrintPath> &out3DPaths);

  /*! \brief Convenience orchestrator for projectPathsTo3D. */
  bool computeProjection(const std::vector<PrintPath> &in2DPaths,
                         std::vector<PrintPath> &out3DPaths) {
    return projectPathsTo3D(in2DPaths, out3DPaths);
  }

  // Color/palette/pipe-mesh helpers and translateMesh/translateGraph/
  // segmentGraph live in zTsUtils, not here. exportToOBJ lives on zTs3DP.
};

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
// clang-format off
#include <source/zToolsets/slicer/zSlicer_Structure.cpp>
#include <source/zToolsets/slicer/zSlicer_Unroll.cpp>
#include <source/zToolsets/slicer/zSlicer_Field.cpp>
// clang-format on
#endif

#endif // ZSPACE_TS_GEOMETRY_SLICER_H
