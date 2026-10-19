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

#ifndef ZSPACE_TS_3DP_ZSLICER_H
#define ZSPACE_TS_3DP_ZSLICER_H

#pragma once

#include <headers/base/zSpace_Toolsets.h>
#include <headers/zToolsets/ts3dp/printsyn/zPrintSyn.h>
#include <headers/zToolsets/ts3dp/utils/Utils.h>

#include <string>

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Collection of toolsets for geometry, analysis, and fabrication.
 *  @{
 */

/** \addtogroup zTs3DP
 *	\brief 3D printing toolset: slicer operators, print-path post-processing.
 *  @{
 */

/*! \class zSlicer
 *	\brief Interface every slicer operator (TopologySlicer, HarmonicSlicer)
 * implements so zTs3DP can run either one: slicing (computeSlicer) and then
 * print-path post-processing (computePostProcess), which calls the zPrintSyn
 * functions directly on the operator's own fields, charts and sections.
 */
class ZSPACE_TOOLSETS zSlicer {
public:
  virtual ~zSlicer() = default;

  /*! \brief The input: a mesh (not owned - must outlive every compute call)
   * plus one VertexTag per mesh FACE, in polygon order (Bottom = print
   * floor, Top = print ceiling, Side = wall - e.g. the OBJ's "usemtl" names
   * read by zTsIO::readMeshFromOBJ()), copied. computeBasicStructure()
   * derives the per-vertex tags from them (zTsUtils::vertexTagsFromFaceTags())
   * and stops if they are inconsistent. inFaceWidthIds: the width group per
   * face (side_<n> -> SlicerSettings::BoundaryWidth[n]; empty = group 0) -
   * the boundary is split where the wall face's group changes. */
  virtual void setInputMesh(zObjMesh &inMesh, const zIntArray &inFaceTags,
                            const zIntArray &inFaceWidthIds = zIntArray()) = 0;

  /*! \brief The shared settings set (copied); each operator reads only its
   * own entries - see SlicerSettings. */
  virtual void setSettings(const SlicerSettings &inSettings) = 0;
  virtual SlicerSettings &getSettings() = 0;
  virtual const SlicerSettings &getSettings() const = 0;

  /*! \brief Optional bracing graph, in the same world coordinates as the
   * input mesh. Not owned. */
  virtual void setBracingGraph(zObjGraph &inGraph) = 0;

  /*! \brief Boundary trim-slot anchor (DoubleBoundary), in the same world
   * coordinates as the input mesh, borrowed (must outlive computation): the slot sits
   * SlicerSettings::BoundarySlotOffset * w either side of the anchor's
   * nearest point on the floor boundary (even layers counter-clockwise, odd
   * clockwise, seen from +Z). An explicit slot point is required before computation. */
  virtual void setSlotAnchor(const zPoint &inAnchor) = 0;
  virtual void clearSlotAnchor() = 0;

  virtual bool isLoaded() const = 0;
  virtual const std::string &getStatus() const = 0;

  /*! \brief The operator's first, input-validating step - cheap enough to
   * run on its own (e.g. to inspect the boundary before configuring). */
  virtual bool computeBasicStructure(bool useParallel = true) = 0;

  /*! \brief Runs the operator's whole pipeline up to and including its 2D
   * SDF path/width/track fields. */
  // false serializes this call's existing OpenMP loops.
  virtual bool computeSlicer(bool useParallel = true) = 0;

  /*! \brief Extracts field paths and maps each layer directly to 3D. */
  virtual bool computePrintPath(bool useParallel = true) = 0;

  /*! \brief Extracts paths, then runs PrintSyn resampling, seam alignment,
   * projection and linking. */
  virtual bool computePostProcess(const PrintSynSetting &settings, bool useParallel = true) = 0;

  /*! \brief Per-layer paths from computePrintPath() or computePostProcess(). */
  virtual const PrintPathArray &getPrintPaths2D() const = 0;
  virtual const PrintPathArray &getPrintPaths3D() const = 0;
  virtual const PrintPath &getResultantPrintPath() const = 0;
};

/** @}*/

/** @}*/

} // namespace zSpace

#endif // ZSPACE_TS_3DP_ZSLICER_H
