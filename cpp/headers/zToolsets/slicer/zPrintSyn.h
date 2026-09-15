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

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

/** \addtogroup zToolsets
 *	\brief Toolsets for geometry, fabrication and computation.
 *  @{
 */

/** \addtogroup zTsGeometry
 *	\brief Geometry toolsets.
 *  @{
 */

/*! \struct PrintSynSetting
 *	\brief Configuration settings for 2D toolpath cleanup and resampling.
 */
struct PrintSynSetting {
  float ResampleLength = 0.15f; // arc-length resampling step (mm)
  float FeatureAngleTol =
      90.0f; // geometric feature corner angle threshold (deg)
};

/*! \class zPrintSyn
 *	\brief Toolpath cleanup, corner-preserving resampling, and width-filtering engine.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zPrintSyn {
protected:
  PrintSynSetting i_Settings;
  std::string m_Status;

public:
  zPrintSyn() { m_Status = "Initialized"; }
  ~zPrintSyn() {}

  void setSettings(const PrintSynSetting &inSettings) {
    i_Settings = inSettings;
  }
  const PrintSynSetting &getSettings() const { return i_Settings; }
  PrintSynSetting &getSettings() { return i_Settings; }

  const std::string &getStatus() const { return m_Status; }

  /*! \brief Runs buildSequencePrintPath -> resamplePrintPath -> filterGlitches -> orientAndWeldLayer on every layer's inner/outer track paths, in parallel per layer. */
  bool computePrintSyn(std::vector<std::vector<PrintPath>> &ioPaths) {
    const int totalLayers = int(ioPaths.size());
    if (totalLayers == 0) {
      m_Status = "failed: No print paths available.";
      zTsUtils::printError("[zPrintSyn] " + m_Status);
      return false;
    }

    const auto wallStart = std::chrono::high_resolution_clock::now();
    int totalPaths = 0;

#pragma omp parallel for schedule(dynamic)
    for (int l = 0; l < totalLayers; l++) {
      for (auto &path : ioPaths[l]) {
        if (buildSequencePrintPath(path) && resamplePrintPath(path)) {
          filterGlitches(path);
        }
      }
      orientAndWeldLayer(ioPaths[l]);
#pragma omp atomic
      totalPaths += int(ioPaths[l].size());
    }

    const auto wallEnd = std::chrono::high_resolution_clock::now();
    const double wallTotalMs =
        std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

    std::ostringstream out;
    out << "Post-processed " << totalPaths << " print paths in " << std::fixed
        << std::setprecision(2) << wallTotalMs
        << " ms (sampleLen=" << i_Settings.ResampleLength << "mm)";
    m_Status = out.str();
    zTsUtils::printSuccess("[zPrintSyn] " + m_Status);
    return true;
  }

  /*! \brief Step 1: directly rebuilds ioPath.graph from its existing vertex order (already cleaned/split/ordered upstream during extraction). */
  static bool buildSequencePrintPath(PrintPath &ioPath,
                                     double weldTolerance = 0.08,
                                     double spurLengthThreshold = 1.0,
                                     std::string *outReason = nullptr);

  /*! \brief Step 2: arc-length resamples ioPath.graph, carrying widths/heights/normals along. */
  bool resamplePrintPath(PrintPath &ioPath);

  /*! \brief Step 3: filters localized width glitches in ioPath.widths. */
  static void filterGlitches(PrintPath &ioPath);

  /*! \brief Step 4: orients one layer's 1-2 ring paths into a single closed loop - a lone path is forced CCW; an inner+outer pair has inner forced CCW and placed first, with outer welded onto its end (starting from inner's end point), closing the loop. Replaces ringPaths with the single merged result. */
  static bool orientAndWeldLayer(std::vector<PrintPath> &ringPaths);

  /*! \brief Cuts a gap of width `width` out of ioPath's closed graph, centered on the point of ioPath nearest to nearPoint. The forward-offset cut point becomes the new start, the backward-offset one the new end. */
  static bool breakPrintPath(PrintPath &ioPath, const zPoint &nearPoint,
                             float width);

  /*! \brief Concatenates every layer's single (post-orientAndWeldLayer) PrintPath end-to-start into one continuous bottom-to-top PrintPath. */
  static bool
  chainLayersIntoOnePath(const std::vector<std::vector<PrintPath>> &layeredPaths,
                        PrintPath &outCombined);
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
