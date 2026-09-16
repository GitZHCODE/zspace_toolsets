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
  // Seed layer only: distance to walk from the extracted start along the
  // inner ring before recording it as the mesh-tracked aligned start point.
  float StartWalkDistance = 5.0f;
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

  /*! \brief Runs buildSequencePrintPath -> resamplePrintPath -> filterGlitches -> orientAndWeldLayer on every layer's inner/outer track paths (in parallel per layer), then alignLayerStartPoints across all layers using sectionMeshes2D. */
  bool computePrintSyn(std::vector<std::vector<PrintPath>> &ioPaths,
                       const zObjMeshArray &sectionMeshes2D) {
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

    alignLayerStartPoints(ioPaths, sectionMeshes2D, i_Settings.StartWalkDistance);

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

  /*! \brief Step 4: orients one layer's ring path CCW in place (extraction already produces at most one ring per layer, already inner-before-outer ordered). */
  static bool orientAndWeldLayer(std::vector<PrintPath> &ringPaths);

  /*! \brief Step 5: aligns every layer's ring start point to the same barycentric point on the underlying mesh (tracks the model through taper/twist), constrained to the inner track. Seeds from the first usable layer by walking walkDistance from its extracted start along the inner segment. */
  static bool
  alignLayerStartPoints(std::vector<std::vector<PrintPath>> &ioPaths,
                       const zObjMeshArray &sectionMeshes2D,
                       double walkDistance);
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
