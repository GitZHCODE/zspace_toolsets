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
#include <headers/zToolsets/ts3dp/utils/Utils.h>
#include <headers/zToolsets/ts3dp/zSlicer/HarmonicSlicer.h>
#include <headers/zToolsets/ts3dp/zSlicer/TopologySlicer.h>
#include <headers/zToolsets/ts3dp/zSlicer/zSlicer.h>
#include <headers/zToolsets/ts3dp/printsyn/zPrintSyn.h>

// Both slicer operators are always compiled (only the selected one is
// instantiated); HarmonicSlicer needs TetGen (include path + tetgen.lib) in
// every build that includes this header.

#include <memory>
#include <string>
#include <vector>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Toolsets for geometry, fabrication and computation.
 *  @{
 */

/** \addtogroup zTs3DP
 *	\brief 3D printing toolset: slicer operators, print-path post-processing.
 *  @{
 */

/*! \enum zTsSlicerMode
 *	\brief Which slicer operator zTs3DP runs.
 */
enum class zTsSlicerMode : int {
  Topology = 0, // TopologySlicer: geodesic V-loop slicing + per-layer 2D unroll
  Harmonic = 1  // HarmonicSlicer: harmonic tet-volume slicing
};

/*! \class zTs3DP
 *	\brief Top-level 3D printing pipeline: one slicer operator (TopologySlicer or
 * HarmonicSlicer, picked by setSlicerMode()); the operator slices
 * (computeSlicer) and post-processes its own result into print paths with the
 * zPrintSyn functions (computePostProcess).
 *
 *  Exactly one operator exists at a time: the constructor creates the
 *  default (Topology) one, and every setSlicerMode() call destroys it and
 *  creates the requested one, bound to the stored input mesh, tags and
 *  bracing. The one SlicerSettings set lives here, so it can be edited in any
 *  mode and survives mode switches; it is copied into the operator at the
 *  start of every compute call.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zTs3DP {
protected:
  //--------------------------
  //---- PROTECTED ATTRIBUTES
  //--------------------------

  zObjMesh i_BlockMesh;
  zObjGraph i_BracingGraph;
  zIntArray i_FaceTags;
  zIntArray i_FaceWidthIds;
  // Boundary trim-slot anchor (world coordinates) - see setSlotAnchor().
  zPoint i_SlotPoint = zPoint(0, 0, 0);
  bool _HasSlotPoint = false;
  zTsSlicerMode _Mode = zTsSlicerMode::Topology;
  SlicerSettings _Settings;
  // The operator for _Mode - never null (created by the constructor and
  // recreated by setSlicerMode()).
  std::unique_ptr<zSlicer> _Slicer;
  PrintSynSetting _PrintSynSettings;
  std::string _Status;

  /*! \brief Replaces _Slicer with a new operator for _Mode, bound to the
   * stored input mesh, tags, bracing and slot anchor. */
  void _createSlicer();


public:
  //--------------------------
  //---- CONSTRUCTOR & DESTRUCTOR
  //--------------------------
  zTs3DP() : _Status("Initialized") { _createSlicer(); }
  ~zTs3DP() = default;
  // Not copyable (the operator holds pointers into this object's input
  // mesh/graph). Movable: the moved operator is re-bound to the target's
  // copies (so `ts = zTs3DP();` resets) and the source gets a fresh one.
  zTs3DP(const zTs3DP &) = delete;
  zTs3DP &operator=(const zTs3DP &) = delete;
  zTs3DP(zTs3DP &&other);
  zTs3DP &operator=(zTs3DP &&other);

  //--------------------------
  //---- SET METHODS
  //--------------------------

  /*! \brief Copies mesh and faceTags; the selected operator is bound to
   * the copies. faceTags: one VertexTag (Bottom/Top/Side - see Utils.h) per
   * mesh face in polygon order, required by both operators (they derive the
   * vertex tags) - e.g. the "usemtl" names read with
   * zTsIO::readMeshFromOBJ(). faceWidthIds: the width group per face
   * (side_<n> -> SlicerSettings::BoundaryWidth[n]; empty = group 0).
   * zTs3DP does no file IO. */
  void setInputMesh(zObjMesh &mesh, const zIntArray &faceTags = zIntArray(),
                    const zIntArray &faceWidthIds = zIntArray());

  /*! \brief Copies graph (world coordinates, same as the input mesh); the
   * selected operator is bound to the copy. */
  void setBracingGraph(zObjGraph &graph);

  /*! \brief Boundary trim-slot anchor (DoubleBoundary), world coordinates
   * like the input mesh - e.g. the first point of a point-only OBJ read with
   * zTsIO::readPointsFromOBJ(). Both operators put the slot
   * SlicerSettings::BoundarySlotOffset * w either side of its nearest point
   * on the floor boundary (even layers counter-clockwise, odd clockwise,
   * seen from +Z). clearSlotAnchor() invalidates the required slot point until set again. */
  void setSlotAnchor(const zPoint &anchor);
  void clearSlotAnchor();
  bool hasSlotAnchor() const { return _HasSlotPoint; }
  const zPoint &getSlotAnchor() const { return i_SlotPoint; }

  /*! \brief Destroys the current operator (and all its results, also when
   * mode is unchanged), creates a new one for mode bound to the stored
   * inputs. */
  void setSlicerMode(zTsSlicerMode mode);

  //--------------------------
  //---- MAIN COMPUTE METHODS
  //--------------------------

  /*! \brief Runs only the selected operator's first, input-validating step.
   * Topology: floor cap, U-loop and V-loops. Harmonic: weld, tag and
   * tetrahedralize. computeSlicer() re-runs it (safe/idempotent on the same
   * inputs). */
  bool computeBasicStructure(bool useParallel = true);

  /*! \brief Runs the selected operator's whole pipeline up to its 2D fields. */
  bool computeSlicer(bool useParallel = true);

  /*! \brief Extracts and directly maps the selected operator's field paths. */
  bool computePrintPath(bool useParallel = true);

  /*! \brief Runs the selected operator's computePostProcess() with
   * getPrintSynSetting(): 2D/3D print paths per layer and the resultant
   * bottom-to-top path (get2DPrintPaths()/get3DPrintPaths()/
   * getResultantPrintPath()). */
  bool computePostProcess(bool useParallel = true);

  //--------------------------
  //---- GET METHODS : pipeline (mode-independent)
  //--------------------------
  inline zTsSlicerMode getSlicerMode() const { return _Mode; }

  /*! \brief The operator for the current mode (replaced by every
   * setSlicerMode() call - don't keep the reference across one). */
  inline zSlicer &getActiveSlicer() { return *_Slicer; }
  inline const zSlicer &getActiveSlicer() const { return *_Slicer; }

  inline bool isLoaded() const { return _Slicer->isLoaded(); }
  inline const std::string &getStatus() const { return _Status; }

  inline PrintSynSetting &getPrintSynSetting() { return _PrintSynSettings; }
  inline const PrintSynSetting &getPrintSynSetting() const {
    return _PrintSynSettings;
  }

  /*! \brief The selected operator's computePostProcess() results (empty
   * until it succeeds; replaced by every setSlicerMode() call). */
  inline const PrintPathArray &get2DPrintPaths() const {
    return _Slicer->getPrintPaths2D();
  }
  inline const PrintPathArray &get3DPrintPaths() const {
    return _Slicer->getPrintPaths3D();
  }
  inline const PrintPath &getResultantPrintPath() const {
    return _Slicer->getResultantPrintPath();
  }

  inline zObjMesh &getInputMesh() { return i_BlockMesh; }
  inline const zObjMesh &getInputMesh() const { return i_BlockMesh; }

  inline const zIntArray &getFaceTags() const { return i_FaceTags; }
  inline const zIntArray &getFaceWidthIds() const { return i_FaceWidthIds; }

  inline zObjGraph &getBracingGraph() { return i_BracingGraph; }
  inline const zObjGraph &getBracingGraph() const { return i_BracingGraph; }

  //--------------------------
  //---- GET METHODS : settings (independent of the operator instance)
  //--------------------------
  /*! \brief The one settings set for both operators (shared entries plus
   * T* Topology-only and H* Harmonic-only ones - see SlicerSettings). */
  inline SlicerSettings &getSettings() { return _Settings; }
  inline const SlicerSettings &getSettings() const { return _Settings; }

  //--------------------------
  //---- GET METHODS : Topology operator (TopologySlicer)
  //--------------------------
  // Null / empty unless Topology is the selected mode.

  TopologySlicer *getTopologySlicer();
  const TopologySlicer *getTopologySlicer() const;

  /*! \brief Number of width groups the input's side faces use (largest
   * side_<n> + 1), i.e. the length SlicerSettings::BoundaryWidth needs - in
   * both modes, straight from the face width ids (no compute needed). */
  int getWidthGroupCount() const;

  zObjMeshArray &getIsoMeshes();
  zObjMeshArray &getIsoMeshes2D();
  zObjMeshScalarField &getScalarField();
  std::vector<zScalarArray> &getPathFields();
  std::vector<zScalarArray> &getWidthFields();
  std::vector<zScalarArray> &getTrackFields();

  //--------------------------
  //---- GET METHODS : Harmonic operator (HarmonicSlicer)
  //--------------------------
  // Null unless Harmonic is the selected mode.

  HarmonicSlicer *getHarmonicSlicer();
  const HarmonicSlicer *getHarmonicSlicer() const;
};

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include <source/zToolsets/ts3dp/zTs3DP.cpp>
#endif

#endif // ZSPACE_TS_GEOMETRY_3DP_H
