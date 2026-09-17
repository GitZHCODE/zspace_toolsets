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
  int StartCornerVertexID = 21;
  float LayerSpacing = 1.0f;
  bool TopoSplit = true;
  bool AngleSplit = true;
  float AngleSplitTol = 90.0f; // degrees (corner angle < threshold)

  // SDF & Toolpath Parameters
  bool DoubleBoundary = false;
  std::vector<float> BoundaryWidth = {0.5f};
  float BracingWidth = 0.5f;
  float OverlapWidth = 0.0f;
  float BracingSlotEven = 0.5f;
  float BracingSlotOdd = 0.5f;
  int BoundarySlotID = 1;
  float BoundarySlotEven = 0.45f;
  float BoundarySlotOdd = 0.55f;

  float FieldScaleUp = 1.05f;
  float FieldCellSize = 0.05f; // Target square grid cell size (default 0.05)
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
  std::vector<int> m_BoundarySplitIndices;
  ParameterGraph m_BracingParams;
  std::string m_Status;

  //---- OUTPUTS
  // sectionMeshes/sectionMeshes2D/pathFields/widthFields/trackFields/the
  // scalar field topology are owned by zTs3DP instead (passed in/out of
  // computeSlice/computeUnroll/computeField as parameters) - they're the
  // data shared with zPrintSyn's extraction/projection steps, so the
  // orchestrator holds them rather than zSlicer.
  zObjGraphArray sectionGraphs;
  zObjGraphArray bracingGraphs;
  zObjGraphArray sectionGraphs2D;
  zObjGraphArray bracingGraphs2D;
  zObjGraphArray trimSlotGraphs2D;

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
  int getLayerCount() const { return int(sectionGraphs.size()); }

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

  zObjGraphArray &getSectionGraphs() { return sectionGraphs; }
  const zObjGraphArray &getSectionGraphs() const { return sectionGraphs; }
  zObjGraphArray &getBracingGraphs() { return bracingGraphs; }
  const zObjGraphArray &getBracingGraphs() const { return bracingGraphs; }

  zObjGraphArray &get2DSectionGraphs() { return sectionGraphs2D; }
  const zObjGraphArray &get2DSectionGraphs() const { return sectionGraphs2D; }
  zObjGraphArray &get2DBracingGraphs() { return bracingGraphs2D; }
  const zObjGraphArray &get2DBracingGraphs() const { return bracingGraphs2D; }

  zObjGraphArray &getTrimSlotGraphs2D() { return trimSlotGraphs2D; }
  const zObjGraphArray &getTrimSlotGraphs2D() const {
    return trimSlotGraphs2D;
  }

  const std::vector<int> &getBoundarySplitIndices() const {
    return m_BoundarySplitIndices;
  }
  const ParameterGraph &getBracingParams() const { return m_BracingParams; }

  //--------------------------
  //---- MAIN COMPUTE METHODS
  //--------------------------

  /*! \brief Corner extraction, V-loops, U-mesh, U-loop, bracing parameterization. */
  bool computeBasicStructure();

  /*! \brief Calculates layer meshes (output via outSectionMeshes, owned by the caller), boundary section graphs, layer bracing graphs. */
  bool computeSlice(zObjMeshArray &outSectionMeshes);

  /*! \brief 2D unrolling of sectionMeshes (input) into outSectionMeshes2D (output, owned by the caller), plus boundary/bracing graphs. */
  bool computeUnroll(const zObjMeshArray &sectionMeshes,
                     zObjMeshArray &outSectionMeshes2D);

  /*! \brief 2D SDF scalar fields only - stops at outScalarField/outPathFields/outWidthFields/outTrackFields, no isocontour extraction (that's zPrintSyn::extractTrackPaths now). sectionMeshes2D is the computeUnroll output, passed back in. */
  bool computeField(const zObjMeshArray &sectionMeshes2D,
                    zObjMeshScalarField &outScalarField,
                    std::vector<zScalarArray> &outPathFields,
                    std::vector<zScalarArray> &outWidthFields,
                    std::vector<zScalarArray> &outTrackFields);

  // Color/palette/pipe-mesh helpers and translateMesh/translateGraph/
  // segmentGraph live in zTsUtils, not here. exportToOBJ lives on zTs3DP.
  // Isocontour extraction and 3D projection live on zPrintSyn now
  // (extractTrackPaths / computeProjection) - zSlicer stops at field
  // computation so the two classes' steps don't interleave.
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
