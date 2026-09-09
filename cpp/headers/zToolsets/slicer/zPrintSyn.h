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
  float sampleLength = 0.15f; // arc-length resampling step (mm)
  float featureAngleThreshold =
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

  /*! \brief Runs buildSequencePrintPath -> resamplePrintPath -> filterGlitches on every path, in parallel. */
  bool computePrintSyn(std::vector<PrintPath> &ioPaths) {
    const int totalLayers = int(ioPaths.size());
    if (totalLayers == 0) {
      m_Status = "failed: No print paths available.";
      zTsUtils::printError("[zPrintSyn] " + m_Status);
      return false;
    }

    const auto wallStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int l = 0; l < totalLayers; l++) {
      if (buildSequencePrintPath(ioPaths[l]) && resamplePrintPath(ioPaths[l])) {
        filterGlitches(ioPaths[l]);
      }
    }

    const auto wallEnd = std::chrono::high_resolution_clock::now();
    const double wallTotalMs =
        std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

    std::ostringstream out;
    out << "Post-processed " << totalLayers << " print paths in " << std::fixed
        << std::setprecision(2) << wallTotalMs
        << " ms (sampleLen=" << i_Settings.sampleLength << "mm)";
    m_Status = out.str();
    zTsUtils::printSuccess("[zPrintSyn] " + m_Status);
    return true;
  }

  /*! \brief Step 1: prunes spurs, welds gaps, rewrites ioPath as the single cleaned ordered chain/loop. */
  static bool buildSequencePrintPath(PrintPath &ioPath,
                                     double weldTolerance = 0.08,
                                     double spurLengthThreshold = 1.0,
                                     std::string *outReason = nullptr);

  /*! \brief Step 2: arc-length resamples ioPath.graph, carrying widths/heights/normals along. */
  bool resamplePrintPath(PrintPath &ioPath);

  /*! \brief Step 3: filters localized width glitches in ioPath.widths. */
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
