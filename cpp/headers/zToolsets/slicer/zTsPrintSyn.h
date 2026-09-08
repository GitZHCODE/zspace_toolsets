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
// zTsSDFSlicer and zTsPrintSyn can both use it without depending on
// each other.

/*! \struct PrintSynSetting
 *	\brief Configuration settings for 2D toolpath cleanup and resampling.
 * Consumed only by zTsPrintSyn.
 */
struct PrintSynSetting {
  float sampleLength = 0.15f; // arc-length resampling step (mm)
  float featureAngleThreshold =
      90.0f; // geometric feature corner angle threshold (deg)
};

/*! \class zTsPrintSyn
 *	\brief Toolpath synthesis, cleanup, corner preservation, uniform resampling, and width filtering engine.
 *	\details Fully decoupled from zTsSDFSlicer: it only depends on zTsUtils.h
 *           (for the shared PrintPath type and shared vector/topology helpers)
 *           and never includes zTsSDFSlicer.h.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zTsPrintSyn {
protected:
  /*!	\brief toolpath cleanup/resampling configuration settings  */
  PrintSynSetting i_Settings;

  std::string m_Status;

public:
  zTsPrintSyn();
  ~zTsPrintSyn();

  /*! \brief Sets toolpath cleanup/resampling configuration settings. */
  void setSettings(const PrintSynSetting &inSettings) {
    i_Settings = inSettings;
  }
  const PrintSynSetting &getSettings() const { return i_Settings; }
  PrintSynSetting &getSettings() { return i_Settings; }

  const std::string &getStatus() const { return m_Status; }

  /*! \brief Processes a single PrintPath (cleans graph, resamples arc-length, filters width glitches). */
  bool processPath(PrintPath &ioPath);

  /*! \brief Processes all PrintPaths in parallel. */
  bool processAllPaths(std::vector<PrintPath> &ioPaths);

  /*! \brief Processes a single layer's 2D contour graph and vertex line widths.
   *
   *	\param[in,out]		ioGraph					- 2D contour graph to clean, sort and resample.
   *	\param[in,out]		ioWeights				- per-vertex extrusion line width weights.
   *	\return				bool					- true if successfully processed.
   */
  bool processLayer(zObjGraph &ioGraph, zFloatArray &ioWeights);

  /*! \brief Prunes spurs, welds gaps, and extracts continuous ordered loop/chain sequence. */
  static bool buildClosedContourSequence(const zObjGraph &graph,
                                         std::vector<int> &sequence,
                                         bool &closed,
                                         std::string *outReason = nullptr,
                                         double weldTolerance = 0.08,
                                         double spurLengthThreshold = 1.0);

  /*! \brief Detects geometric corner feature points preserving sharp angles. */
  static void detectFeatureAngles(const zPointArray &positions,
                                  const std::vector<int> &sequence,
                                  float angleThresholdDeg,
                                  std::vector<int> &outFeatureFlags);

  /*! \brief Performs uniform arc-length resampling preserving detected feature corners. */
  static void resampleContourWithWeights(const zPointArray &positions,
                                        const zFloatArray &weights,
                                        const std::vector<int> &sequence,
                                        const std::vector<int> &featureFlags,
                                        bool closed, float sampleLength,
                                        zPointArray &outPoints,
                                        zFloatArray &outWeights,
                                        zIntArray &outEdges);

  /*! \brief Filters out localized 1-vertex extrusion width glitches and noise. */
  static void filterWidthGlitches(zFloatArray &widths, bool closed);
};

/** @}*/

/** @}*/

} // namespace zSpace

#if defined(ZSPACE_TOOLSETS_STATIC_LIBRARY) ||                                 \
    defined(ZSPACE_TOOLSETS_DYNAMIC_LIBRARY)
// All defined OK so do nothing
#else
#include <source/zToolsets/slicer/zTsPrintSyn.cpp>
#endif

#endif // ZSPACE_TS_GEOMETRY_PRINTSYN_H
