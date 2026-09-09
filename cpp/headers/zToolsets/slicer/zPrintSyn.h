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

#ifndef ZSPACE_TS_GEOMETRY_PRINTSYN_H
#define ZSPACE_TS_GEOMETRY_PRINTSYN_H

#pragma once

#include <headers/base/zSpace_Toolsets.h>
#include <headers/zToolsets/slicer/zTsUtils.h>
#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/zInterface/functionsets/zFnGraph.h>
#include <zspace/zInterface/functionsets/zFnMesh.h>
#include <zspace/zInterface/objects/zObjGraph.h>
#include <zspace/zInterface/objects/zObjMesh.h>

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

// NOTE: PrintPath is a shared type declared in zTsUtils.h so that
// zSlicer and zPrintSyn can both use it without depending on
// each other.

/*! \struct PrintSynSetting
 *	\brief Configuration settings for 2D toolpath cleanup and resampling.
 * Consumed only by zPrintSyn.
 */
struct PrintSynSetting {
  float sampleLength = 0.15f; // arc-length resampling step (mm)
  float featureAngleThreshold =
      90.0f; // geometric feature corner angle threshold (deg)
};

/*! \class zPrintSyn
 *	\brief Toolpath synthesis, cleanup, corner preservation, uniform resampling, and width filtering engine.
 *	\details Fully decoupled from zSlicer: it only depends on zTsUtils.h
 *           (for the shared PrintPath type and shared vector/topology helpers)
 *           and never includes zSlicer.h.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zPrintSyn {
protected:
  /*!	\brief toolpath cleanup/resampling configuration settings  */
  PrintSynSetting i_Settings;

  std::string m_Status;

public:
  /*! \brief Default constructor. */
  zPrintSyn() { m_Status = "Initialized"; }

  /*! \brief Default destructor. */
  ~zPrintSyn() {}

  /*! \brief Sets toolpath cleanup/resampling configuration settings. */
  void setSettings(const PrintSynSetting &inSettings) {
    i_Settings = inSettings;
  }
  const PrintSynSetting &getSettings() const { return i_Settings; }
  PrintSynSetting &getSettings() { return i_Settings; }

  const std::string &getStatus() const { return m_Status; }

  /*! \brief Processes all PrintPaths in parallel: runs
   * buildSequencePrintPath -> resamplePrintPath -> filterGlitches on each,
   * in sequence. Each step takes/returns the same PrintPath, carrying
   * .widths/.heights/.normals along with .graph throughout. Works for both
   * 2D paths (.heights/.normals empty) and 3D paths (all four populated),
   * since empty channels are simply skipped.
   */
  bool processAllPaths(std::vector<PrintPath> &ioPaths);

  /*! \brief Step 1: prunes spurs, welds gaps, extracts the longest
   * continuous loop/chain, and rewrites ioPath.graph as that single clean,
   * ordered sequence (widths/heights/normals reordered to match). Independent
   * step - takes and returns a complete PrintPath.
   *	\param[in,out]		ioPath					- print path to clean.
   *	\return				bool					- true if a valid contour was extracted.
   */
  static bool buildSequencePrintPath(PrintPath &ioPath,
                                     double weldTolerance = 0.08,
                                     double spurLengthThreshold = 1.0,
                                     std::string *outReason = nullptr);

  /*! \brief Step 2: performs uniform arc-length resampling preserving
   * detected feature corners, carrying widths/heights/normals along with
   * positions (heights/normals skipped when empty). Assumes ioPath.graph is
   * already a single ordered chain/loop (e.g. from
   * buildSequencePrintPath). Reads sampleLength/featureAngleThreshold
   * from this instance's settings. Independent step - takes and returns a
   * complete PrintPath.
   *	\param[in,out]		ioPath					- print path to resample.
   *	\return				bool					- true if successfully resampled.
   */
  bool resamplePrintPath(PrintPath &ioPath);

  /*! \brief Step 3: filters out localized 1-vertex extrusion width glitches
   * and noise in ioPath.widths (closed/open inferred from ioPath.graph).
   * Independent step - takes and returns a complete PrintPath.
   *	\param[in,out]		ioPath					- print path whose widths to filter.
   */
  static void filterGlitches(PrintPath &ioPath);
};

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include <source/zToolsets/slicer/zPrintSyn.cpp>
#endif

#endif // ZSPACE_TS_GEOMETRY_PRINTSYN_H
