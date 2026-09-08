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

#include <sstream>

namespace zSpace {

//---- MAIN COMPUTE METHODS: STEP 4 (FINALIZE 2D PRINT PATHS)

ZSPACE_INLINE bool zTsSDFSlicer::computePath() {
  // computeField() already produced _2DPrintPaths (isocontour graph +
  // per-vertex widths per layer). This step only validates/exposes that
  // result - projecting to 3D is a separate, explicit step
  // (projectPathsTo3D, in zTsSDFSlicer_Unroll.cpp) so that 2D
  // post-processing can run on _2DPrintPaths in between without this step
  // redoing a 3D projection on data that hasn't changed.
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

} // namespace zSpace
