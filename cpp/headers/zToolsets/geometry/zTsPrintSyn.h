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

/*! \struct PrintPath
 *	\brief Encapsulates a complete printing toolpath including its graph, line widths, layer heights, and normals.
 */
struct PrintPath {
  zObjGraph graph;
  zFloatArray widths;
  zFloatArray heights;
  zVectorArray normals;

  PrintPath() = default;
  PrintPath(const zObjGraph &inGraph, const zFloatArray &inWidths = {},
            const zFloatArray &inHeights = {},
            const zVectorArray &inNormals = {})
      : graph(inGraph), widths(inWidths), heights(inHeights),
        normals(inNormals) {}
};

/*! \class zTsPrintSyn
 *	\brief Toolpath synthesis, cleanup, corner preservation, uniform resampling, and width filtering engine.
 *	\details Decoupled post-processing engine that transforms raw isocontours into continuous,
 *           high-quality robotic 3D printing toolpaths.
 *	\since 1.0.0
 */
class ZSPACE_TOOLSETS zTsPrintSyn {
protected:
  std::string m_Status;

public:
  zTsPrintSyn();
  ~zTsPrintSyn();

  const std::string &getStatus() const { return m_Status; }

  /*! \brief Processes a single PrintPath (cleans graph, resamples arc-length, filters width glitches). */
  bool processPath(PrintPath &ioPath, float sampleLength = 0.005f,
                   float featureAngleThreshold = 120.0f);

  /*! \brief Processes all PrintPaths in parallel. */
  bool processAllPaths(std::vector<PrintPath> &ioPaths,
                       float sampleLength = 0.005f,
                       float featureAngleThreshold = 120.0f);

  /*! \brief Processes a single layer's 2D contour graph and vertex line widths.
   *
   *	\param[in,out]		ioGraph					- 2D contour graph to clean, sort and resample.
   *	\param[in,out]		ioWeights				- per-vertex extrusion line width weights.
   *	\param[in]			sampleLength			- arc-length resampling step in mm (default 0.005).
   *	\param[in]			featureAngleThreshold	- feature corner angle threshold in deg (default 120.0).
   *	\return				bool					- true if successfully processed.
   */
  bool processLayer(zObjGraph &ioGraph, zFloatArray &ioWeights,
                    float sampleLength = 0.005f,
                    float featureAngleThreshold = 120.0f);

  /*! \brief Processes all layers' 2D contour graphs in parallel.
   *
   *	\param[in,out]		ioGraphs				- array of 2D contour graphs.
   *	\param[in,out]		ioWeights				- array of per-layer vertex line width weights.
   *	\param[in]			sampleLength			- arc-length resampling step in mm (default 0.005).
   *	\param[in]			featureAngleThreshold	- feature corner angle threshold in deg (default 120.0).
   *	\return				bool					- true if successfully processed.
   */
  bool processAllLayers(std::vector<zObjGraph> &ioGraphs,
                        std::vector<zFloatArray> &ioWeights,
                        float sampleLength = 0.005f,
                        float featureAngleThreshold = 120.0f);

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
#include <source/zToolsets/geometry/zTsPrintSyn.cpp>
#endif

#endif // ZSPACE_TS_GEOMETRY_PRINTSYN_H
