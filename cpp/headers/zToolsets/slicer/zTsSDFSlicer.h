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

#ifndef ZSPACE_TS_GEOMETRY_SDFSLICER_H
#define ZSPACE_TS_GEOMETRY_SDFSLICER_H

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

/*! \class zTsSDFSlicer
 *	\brief A toolset for SDF based geodesic slicing and bracing generation.
 *	\details Fully decoupled from zTsPrintSyn: it only produces/consumes
 *           PrintPath and ParameterGraph data (declared in zTsUtils.h),
 *           and never includes zTsPrintSyn.h.
 */
class ZSPACE_TOOLSETS zTsSDFSlicer {
protected:
  //--------------------------
  //---- PROTECTED ATTRIBUTES
  //--------------------------

  //---- INPUTS (i_)
  /*!	\brief pointer to input mesh Object  */
  zObjMesh *i_InputMesh;

  /*!	\brief pointer to input bracing graph Object  */
  zObjGraph *i_BracingGraph;

  /*!	\brief slicer configuration settings  */
  SlicerSettings i_Settings;

  /*!	\brief container of U-mesh (bottom basic mesh) Object  */
  zObjMesh m_UMesh;

  /*!	\brief container of U-loop boundary graph Object  */
  zObjGraph m_ULoop;

  /*!	\brief container of V longitude graph loops  */
  zObjGraphArray m_VLoops;

  /*!	\brief container of V longitude scalar values per loop  */
  std::vector<zScalarArray> m_VScalar;

  /*! \brief intermediate reusable 2D SDF and scalar field operator  */
  zObjMeshScalarField m_ScalarField;

  /*!	\brief boundary split indices for corner segments  */
  std::vector<int> m_BoundarySplitIndices;

  /*!	\brief mesh barycentric parameterized bracing parameters (2D/3D
   * mesh-based)  */
  ParameterGraph m_BracingParams;

  /*!	\brief status string message  */
  std::string m_Status;

  //---- OUTPUTS
  /*!	\brief container of section mesh objects  */
  zObjMeshArray sectionMeshes;

  /*!	\brief container of section graph objects  */
  zObjGraphArray sectionGraphs;

  /*!	\brief container of layer bracing graph objects  */
  zObjGraphArray bracingGraphs;

  /*!	\brief container of 2D unrolled section mesh objects  */
  zObjMeshArray sectionMeshes2D;

  /*!	\brief container of 2D unrolled section graph objects  */
  zObjGraphArray sectionGraphs2D;

  /*!	\brief container of 2D unrolled layer bracing graph objects  */
  zObjGraphArray bracingGraphs2D;

  /*!	\brief container of 2D SDF path scalar values per layer  */
  std::vector<zScalarArray> pathFields;

  /*!	\brief container of 2D width scalar values per layer  */
  std::vector<zScalarArray> widthFields;

  /*!	\brief container of 2D trim slot graphs  */
  zObjGraphArray trimSlotGraphs;

public:
  //--------------------------
  //---- CONSTRUCTOR & DESTRUCTOR
  //--------------------------

  /*! \brief Default constructor. */
  zTsSDFSlicer() {
    i_InputMesh = nullptr;
    i_BracingGraph = nullptr;
    m_Status = "Ready";
  }

  /*! \brief Overloaded constructor with settings. */
  explicit zTsSDFSlicer(const SlicerSettings &inSettings) : zTsSDFSlicer() {
    i_Settings = inSettings;
  }

  /*! \brief Default destructor. */
  ~zTsSDFSlicer() {}

  //--------------------------
  //---- SET METHODS
  //--------------------------

  /*! \brief Sets slicer configuration settings. */
  void setSettings(const SlicerSettings &inSettings) {
    i_Settings = inSettings;
  }
  const SlicerSettings &getSettings() const { return i_Settings; }
  SlicerSettings &getSettings() { return i_Settings; }

  /*! \brief Sets input mesh object. */
  void setInputMesh(zObjMesh &inMesh) { i_InputMesh = &inMesh; }

  /*! \brief Sets input bracing graph object. */
  void setBracingGraph(zObjGraph &inGraph) { i_BracingGraph = &inGraph; }

  /*! \brief Sets current status message. */
  void setStatus(const std::string &inStatus) { m_Status = inStatus; }
  const std::string &getStatus() const { return m_Status; }

  //--------------------------
  //---- GET METHODS
  //--------------------------

  /*! \brief Returns true if slice mesh is loaded. */
  bool isLoaded() const { return i_InputMesh != nullptr; }

  /*! \brief Returns total computed layer count. */
  int getLayerCount() const { return int(sectionMeshes.size()); }

  /*! \brief Returns input mesh object. */
  zObjMesh &getInputMesh() { return *i_InputMesh; }
  const zObjMesh &getInputMesh() const { return *i_InputMesh; }

  /*! \brief Returns bottom U-mesh and U-loop graph. */
  zObjMesh &getUMesh() { return m_UMesh; }
  const zObjMesh &getUMesh() const { return m_UMesh; }
  zObjGraph &getULoop() { return m_ULoop; }
  const zObjGraph &getULoop() const { return m_ULoop; }

  /*! \brief Returns V-loops (longitude graph loops). */
  zObjGraphArray &getVLoops() { return m_VLoops; }
  const zObjGraphArray &getVLoops() const { return m_VLoops; }

  /*! \brief Returns V longitude scalar values per loop. */
  std::vector<zScalarArray> &getVScalar() { return m_VScalar; }
  const std::vector<zScalarArray> &getVScalar() const { return m_VScalar; }

  /*! \brief Returns section mesh objects. */
  zObjMeshArray &getSectionMeshes() { return sectionMeshes; }
  const zObjMeshArray &getSectionMeshes() const { return sectionMeshes; }
  zObjGraphArray &getSectionGraphs() { return sectionGraphs; }
  const zObjGraphArray &getSectionGraphs() const { return sectionGraphs; }
  zObjGraphArray &getBracingGraphs() { return bracingGraphs; }
  const zObjGraphArray &getBracingGraphs() const { return bracingGraphs; }

  /*! \brief Returns 2D unrolled section mesh & graph objects. */
  zObjMeshArray &get2DSectionMeshes() { return sectionMeshes2D; }
  const zObjMeshArray &get2DSectionMeshes() const { return sectionMeshes2D; }
  zObjGraphArray &get2DSectionGraphs() { return sectionGraphs2D; }
  const zObjGraphArray &get2DSectionGraphs() const { return sectionGraphs2D; }
  zObjGraphArray &get2DBracingGraphs() { return bracingGraphs2D; }
  const zObjGraphArray &get2DBracingGraphs() const { return bracingGraphs2D; }

  /*! \brief Returns 2D SDF path and width scalar values per layer. */
  std::vector<zScalarArray> &getPathFields() { return pathFields; }
  const std::vector<zScalarArray> &getPathFields() const { return pathFields; }
  std::vector<zScalarArray> &getWidthFields() { return widthFields; }
  const std::vector<zScalarArray> &getWidthFields() const {
    return widthFields;
  }

  /*! \brief Returns intermediate reusable 2D scalar field operator. */
  zObjMeshScalarField &getScalarField() { return m_ScalarField; }
  const zObjMeshScalarField &getScalarField() const { return m_ScalarField; }

  /*! \brief Returns 2D trim slot graphs. */
  zObjGraphArray &getTrimSlotGraphs() { return trimSlotGraphs; }
  const zObjGraphArray &getTrimSlotGraphs() const { return trimSlotGraphs; }

  /*! \brief Returns corner split indices and bracing barycentric params. */
  const std::vector<int> &getBoundarySplitIndices() const {
    return m_BoundarySplitIndices;
  }
  const ParameterGraph &getBracingParams() const { return m_BracingParams; }

  //--------------------------
  //---- MAIN COMPUTE METHODS
  //--------------------------

  /*! \brief Computes basic structure: corner 3-halfedge extraction, V-loops,
   * U-mesh, U-loop, and bracing parameterization.
   *
   *	\return				bool
   *	- true if successful.
   */
  bool computeBasicStructure();

  /*! \brief Computes slices: calculates layer meshes, boundary section graphs,
   * and layer bracing graphs from U-mesh.
   *
   *	\return				bool
   *	- true if successful.
   */
  bool computeSlice();

  /*! \brief Computes 2D unrolling for all 3D section meshes, boundary graphs,
   * and bracing graphs.
   *
   *	\return				bool
   *	- true if successful.
   */
  bool computeUnroll();

  /*! \brief Computes 2D SDF scalar fields (path field, width field, and trim
   * slots) for all unrolled layers, and extracts the resulting flat
   * isocontour print paths into outPrintPaths2D.
   *
   *	\param[out]			outPrintPaths2D	- 2D flat isocontour
   * print paths (graph + per-vertex widths), one per layer.
   *	\return				bool
   *	- true if successful.
   */
  bool computeField(std::vector<PrintPath> &outPrintPaths2D);

  /*! \brief Projects 2D PrintPaths back to 3D section meshes, populating 3D
   * graphs, heights, and normals. */
  bool projectPathsTo3D(const std::vector<PrintPath> &in2DPaths,
                        std::vector<PrintPath> &out3DPaths);

  // NOTE: color/palette helpers (colorPathField, colorWidthField, hslToRgb,
  // generateColorPalette), draw/export helpers (drawPathMesh, exportToOBJ)
  // are not members of zTsSDFSlicer - drawing/exporting reads print-path
  // data that only the caller (zTs3DP) owns, so those live at the zTs3DP
  // level instead. Call zTsUtils::... (declared in zTsUtils.h) for the
  // generic color/palette/OBJ-export utilities.

  // NOTE: translateMesh, translateGraph, and segmentGraph are not members of
  // zTsSDFSlicer - they are generic geometry utilities that don't depend on
  // any slicer state. Call them via zTsUtils::... (declared in zTsUtils.h)
  // instead.
};

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
// clang-format off
#include <source/zToolsets/slicer/zTsSDFSlicer_Structure.cpp>
#include <source/zToolsets/slicer/zTsSDFSlicer_Unroll.cpp>
#include <source/zToolsets/slicer/zTsSDFSlicer_Field.cpp>
// clang-format on
#endif

#endif // ZSPACE_TS_GEOMETRY_SDFSLICER_H
