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

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/slicer/zTsSDFSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

namespace zSpace {

using namespace zTsUtils;

//---- MAIN COMPUTE METHODS: STEP 4 (FINALIZE 2D PRINT PATHS)

ZSPACE_INLINE bool zTsSDFSlicer::computePath() {
  // computeField() already produced _2DPrintPaths (isocontour graph +
  // per-vertex widths per layer). This step only validates/exposes that
  // result - projecting to 3D is a separate, explicit step
  // (project2DPrintPathsTo3D) so that 2D post-processing can run on
  // _2DPrintPaths in between without this step redoing a 3D projection on
  // data that hasn't changed.
  if (_2DPrintPaths.empty()) {
    m_Status = "failed: No 2D print paths available. Call computeField() first.";
    printError("[computePath] " + m_Status);
    return false;
  }

  std::ostringstream out;
  out << "2D print paths ready: " << _2DPrintPaths.size() << " layers.";
  m_Status = out.str();
  printSuccess("[computePath] " + m_Status);
  return true;
}

//---- MAIN COMPUTE METHODS: STEP 5 (PROJECT 2D PRINT PATHS TO 3D)

ZSPACE_INLINE bool zTsSDFSlicer::project2DPrintPathsTo3D(
    const std::vector<PrintPath> &in2DPaths,
    std::vector<PrintPath> &out3DPaths) {
  const int totalLayers = int(in2DPaths.size());
  if (totalLayers == 0 || _2DSectionMeshes.empty() || sectionMeshes.empty()) {
    m_Status = "failed: No 2D print paths or section meshes available for "
               "projection.";
    printError("[project2DPrintPathsTo3D] " + m_Status);
    return false;
  }

  out3DPaths.assign(totalLayers, PrintPath());
  _3DPrintPaths.assign(totalLayers, PrintPath());

  zFloatArray meshVertexHeights;
  if (!m_VScalar.empty() && totalLayers > 1) {
    meshVertexHeights.resize(m_VScalar.size());
    const float denom = float(totalLayers - 1);
    for (size_t k = 0; k < m_VScalar.size(); k++) {
      meshVertexHeights[k] = (!m_VScalar[k].empty())
                                 ? float(m_VScalar[k].back() / denom)
                                 : i_Settings.layerSpacing;
    }
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    out3DPaths[l].graph = in2DPaths[l].graph;
    out3DPaths[l].widths = in2DPaths[l].widths;

    zFnGraph fnContour(out3DPaths[l].graph);
    if (fnContour.numVertices() > 0 && l < int(_2DSectionMeshes.size()) &&
        l < int(sectionMeshes.size())) {
      barycentricProjectionTriMesh(
          out3DPaths[l].graph, _2DSectionMeshes[l], sectionMeshes[l],
          &out3DPaths[l].normals, &meshVertexHeights, &out3DPaths[l].heights);
    }

    _3DPrintPaths[l] = out3DPaths[l];
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Projected " << out3DPaths.size() << " 3D PrintPaths in "
      << std::fixed << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  printSuccess("[project2DPrintPathsTo3D] " + m_Status);
  return true;
}

} // namespace zSpace
