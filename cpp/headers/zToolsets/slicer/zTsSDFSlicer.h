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

#include <iostream>
#include <mutex>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

// NOTE: SlicerSettings, MeshBarycentricGraph and PrintPath are shared types
// declared in zTsUtils.h so that zTsSDFSlicer and zTsPrintSyn can both use
// them without depending on each other.

/*! \class zTsSDFSlicer
 *	\brief A toolset for SDF based geodesic slicing and bracing generation.
 *	\details Fully decoupled from zTsPrintSyn: it only produces/consumes
 *           PrintPath and MeshBarycentricGraph data (declared in zTsUtils.h),
 *           and never includes zTsPrintSyn.h.
 */
class ZSPACE_TOOLSETS zTsSDFSlicer {
private:
  //--------------------------
  //---- PRIVATE LOGGING HELPERS
  //--------------------------

  static inline std::mutex &_getConsoleLogMutex() {
    static std::mutex mtx;
    return mtx;
  }

  static inline void printSuccess(const std::string &msg) {
    std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
    std::cout << "\033[1;32m" << msg << "\033[0m\n" << std::flush;
  }

  static inline void printError(const std::string &msg) {
    std::lock_guard<std::mutex> lock(_getConsoleLogMutex());
    std::cerr << "\033[1;31m" << msg << "\033[0m\n" << std::flush;
  }

protected:
  //--------------------------
  //---- PROTECTED ATTRIBUTES
  //--------------------------

  /*!	\brief core utilities Object  */
  zUtilsCore coreUtils;

  //---- INPUTS (i_)
  /*!	\brief pointer to input mesh Object  */
  zObjMesh *i_InputMesh;

  /*!	\brief pointer to input bracing graph Object  */
  zObjGraph *i_BracingGraph;

  /*!	\brief slicer configuration settings  */
  SlicerSettings i_Settings;

  //---- MIDDLE / REUSABLE DATA (m_)
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
  MeshBarycentricGraph m_BracingParams;

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
  zObjMeshArray _2DSectionMeshes;

  /*!	\brief container of 2D unrolled section graph objects  */
  zObjGraphArray _2DSectionGraphs;

  /*!	\brief container of 2D unrolled layer bracing graph objects  */
  zObjGraphArray _2DBracingGraphs;

  /*!	\brief container of 2D SDF path scalar values per layer  */
  std::vector<zScalarArray> pathFields;

  /*!	\brief container of 2D width scalar values per layer  */
  std::vector<zScalarArray> widthFields;

  /*!	\brief container of 2D trim slot graphs  */
  zObjGraphArray trimSlotGraphs;

  /*!	\brief container of 2D flat isocontour toolpath graphs  */
  zObjGraphArray _2DContourGraphs;

  /*!	\brief container of 2D/3D isocontour vertex widths per layer  */
  std::vector<zFloatArray> contourWeights;

  /*!	\brief container of 3D contour vertex layer heights per layer  */
  std::vector<zFloatArray> contourHeights;

  /*!	\brief container of 3D contour vertex normal vectors per layer  */
  std::vector<zVectorArray> contourNormals;

  /*!	\brief container of 3D projected contour/toolpath graphs  */
  zObjGraphArray contourGraphs;

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

  /*! \brief Sets target square cell size for SDF field meshes. */
  void setFieldCellSize(float cellSize) {
    i_Settings.fieldCellSize = (cellSize > 1e-6f) ? cellSize : 0.05f;
  }
  float getFieldCellSize() const { return i_Settings.fieldCellSize; }

  /*! \brief Sets input mesh object. */
  void setInputMesh(zObjMesh &inMesh) { i_InputMesh = &inMesh; }
  inline void setSliceMesh(zObjMesh &inMesh) { setInputMesh(inMesh); }

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
  inline int getNumSections() const { return getLayerCount(); }

  /*! \brief Returns input mesh object. */
  zObjMesh &getInputMesh() { return *i_InputMesh; }
  const zObjMesh &getInputMesh() const { return *i_InputMesh; }
  inline zObjMesh &getSliceMesh() { return getInputMesh(); }
  inline const zObjMesh &getSliceMesh() const { return getInputMesh(); }

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
  zObjMeshArray &get2DSectionMeshes() { return _2DSectionMeshes; }
  const zObjMeshArray &get2DSectionMeshes() const { return _2DSectionMeshes; }
  zObjGraphArray &get2DSectionGraphs() { return _2DSectionGraphs; }
  const zObjGraphArray &get2DSectionGraphs() const { return _2DSectionGraphs; }
  zObjGraphArray &get2DBracingGraphs() { return _2DBracingGraphs; }
  const zObjGraphArray &get2DBracingGraphs() const { return _2DBracingGraphs; }

  /*! \brief Returns 2D SDF path and width scalar values per layer. */
  std::vector<zScalarArray> &getPathFields() { return pathFields; }
  const std::vector<zScalarArray> &getPathFields() const {
    return pathFields;
  }
  std::vector<zScalarArray> &getWidthFields() { return widthFields; }
  const std::vector<zScalarArray> &getWidthFields() const {
    return widthFields;
  }

  /*! \brief Returns intermediate reusable 2D scalar field operator. */
  zObjMeshScalarField &getScalarField() { return m_ScalarField; }
  const zObjMeshScalarField &getScalarField() const { return m_ScalarField; }
  inline zObjMeshScalarField &getFieldOperator() { return m_ScalarField; }
  inline const zObjMeshScalarField &getFieldOperator() const {
    return m_ScalarField;
  }
  inline zObjMeshScalarField &getPathFieldOperator() { return m_ScalarField; }
  inline const zObjMeshScalarField &getPathFieldOperator() const {
    return m_ScalarField;
  }
  inline zObjMeshScalarField &getWidthFieldOperator() { return m_ScalarField; }
  inline const zObjMeshScalarField &getWidthFieldOperator() const {
    return m_ScalarField;
  }

  /*! \brief Returns 2D trim slot and isocontour toolpath graphs. */
  zObjGraphArray &getTrimSlotGraphs() { return trimSlotGraphs; }
  const zObjGraphArray &getTrimSlotGraphs() const { return trimSlotGraphs; }
  zObjGraphArray &get2DContourGraphs() { return _2DContourGraphs; }
  const zObjGraphArray &get2DContourGraphs() const { return _2DContourGraphs; }

  /*! \brief Returns contour vertex line widths per layer. */
  inline std::vector<zFloatArray> &get2DContourWeights() {
    return getContourWeights();
  }
  inline const std::vector<zFloatArray> &get2DContourWeights() const {
    return getContourWeights();
  }
  std::vector<zFloatArray> &getContourWeights() { return contourWeights; }
  const std::vector<zFloatArray> &getContourWeights() const {
    return contourWeights;
  }

  /*! \brief Returns 3D contour heights, normals and projected toolpath graphs.
   */
  std::vector<zFloatArray> &getContourHeights() { return contourHeights; }
  const std::vector<zFloatArray> &getContourHeights() const {
    return contourHeights;
  }
  std::vector<zVectorArray> &getContourNormals() { return contourNormals; }
  const std::vector<zVectorArray> &getContourNormals() const {
    return contourNormals;
  }
  zObjGraphArray &getContourGraphs() { return contourGraphs; }
  const zObjGraphArray &getContourGraphs() const { return contourGraphs; }

  /*! \brief Returns corner split indices and bracing barycentric params. */
  const std::vector<int> &getBoundarySplitIndices() const {
    return m_BoundarySplitIndices;
  }
  const MeshBarycentricGraph &getBracingParams() const {
    return m_BracingParams;
  }

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
   * slots) for all unrolled layers.
   *
   *	\return				bool
   *	- true if successful.
   */
  bool computeField();

  /*! \brief Computes 2D isocontour toolpaths and projects them back to 3D
   * section meshes.
   *
   *	\return				bool
   *	- true if successful.
   */
  bool computePath();

  /*! \brief Extracts raw 2D PrintPaths from 2D isocontour graphs and width weights. */
  bool extract2DPrintPaths(std::vector<PrintPath> &outPaths) const;

  /*! \brief Projects 2D PrintPaths back to 3D section meshes, populating 3D graphs, heights, and normals. */
  bool project2DPrintPathsTo3D(const std::vector<PrintPath> &in2DPaths,
                               std::vector<PrintPath> &out3DPaths);

  /*! \brief Exports each layer's section mesh and 3D toolpath graph to
   * directory in OBJ format.
   *
   *	\param[in]			directoryPath			- export
   * folder path (defaults to Desktop/SDF_Slicer_Export if empty).
   *	\return				bool
   *	- true if successful.
   */
  bool exportToOBJ(const std::string &directoryPath = "");

  /*! \brief Computes vertex colors for a scalar field mesh based on path SDF scalar values.
   *
   *	\param[in,out]		mesh		- target mesh whose vertex colors are updated.
   *	\param[in]			values		- scalar values per vertex.
   */
  static void colorPathField(zObjMesh &mesh, const zScalarArray &values);

  /*! \brief Computes vertex color array based on path SDF scalar values.
   *
   *	\param[in]			values		- scalar values per vertex.
   *	\param[out]			outColors	- output vertex colors.
   */
  static void colorPathField(const zScalarArray &values, std::vector<zColor> &outColors);

  /*! \brief Computes discrete vertex colors for a width field mesh grouped by unique path width values.
   *
   *	\param[in,out]		mesh		- target mesh whose vertex colors are updated.
   *	\param[in]			values		- width scalar values per vertex.
   */
  static void colorWidthField(zObjMesh &mesh, const zScalarArray &values);

  /*! \brief Converts HSL color values (H in degrees [0, 360), S in [0, 1], L in [0, 1]) to zColor RGB. */
  static zColor hslToRgb(float h, float s, float l);

  /*! \brief Generates a distinct shuffled HSL color palette for K categories.
   *
   *	\param[in]			numColors		- number of distinct categories (K).
   *	\param[out]			outColors		- generated color palette.
   *	\param[in]			seed			- random seed for shuffling (default 1337).
   *	\param[in]			saturation		- HSL saturation (default 0.85).
   *	\param[in]			lightness		- HSL lightness (default 0.50).
   */
  static void generateColorPalette(int numColors, std::vector<zColor> &outColors,
                                   unsigned int seed = 1337,
                                   float saturation = 0.85f,
                                   float lightness = 0.50f);

  /*! \brief Computes discrete vertex color array grouped by unique path width values.
   *
   *	\param[in]			values		- width scalar values per vertex.
   *	\param[out]			outColors	- output vertex colors.
   */
  static void colorWidthField(const zScalarArray &values, std::vector<zColor> &outColors);

  /*! \brief Builds 3D extruded rectangular pipe mesh for a specific toolpath layer.
   *
   *	\param[in]			layerIdx		- layer index.
   *	\param[out]			outMesh			- generated extruded pipe mesh.
   */
  void drawPathMesh(int layerIdx, zObjMesh &outMesh) const;

  /*! \brief Generic 3D translation for all vertex positions of a mesh. */
  static void translateMesh(zObjMesh &mesh, const zVector &offset);

  /*! \brief Generic 3D translation for all vertex positions of a graph. */
  static void translateGraph(zObjGraph &graph, const zVector &offset);

  /*! \brief Segments a graph at specified split vertex indices. If no indices are provided, outputs original graph. */
  static bool segmentGraph(const zObjGraph &sourceBoundaryGraph,
                           const std::vector<int> &splitIndices,
                           std::vector<zObjGraph> &outBoundarySegments);
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
#include <source/zToolsets/slicer/zTsSDFSlicer_Path.cpp>
#include <source/zToolsets/slicer/zTsSDFSlicer_Export.cpp>
// clang-format on
#endif

#endif // ZSPACE_TS_GEOMETRY_SDFSLICER_H
