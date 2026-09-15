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

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <headers/zToolsets/slicer/zPrintSyn.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace zSpace {

using namespace zTsUtils;

//--------------------------
//---- PRIVATE HELPERS
//--------------------------

namespace {

// Flags sequence indices whose turn angle is below the feature threshold.
void _detectFeatureAngles(const zPointArray &positions,
                          const std::vector<int> &sequence,
                          float angleThresholdDeg,
                          std::vector<int> &outFeatureFlags) {
  const int N = static_cast<int>(sequence.size());
  outFeatureFlags.assign(N, 0);
  if (N < 3)
    return;

  const double thresholdRad =
      static_cast<double>(angleThresholdDeg) * (3.14159265358979323846 / 180.0);

  for (int i = 0; i < N; i++) {
    int prevId = sequence[(i - 1 + N) % N];
    int curId = sequence[i];
    int nextId = sequence[(i + 1) % N];

    zVector v0 = vecSub(positions[prevId], positions[curId]);
    zVector v1 = vecSub(positions[nextId], positions[curId]);
    const double len0 = v0.length();
    const double len1 = v1.length();
    if (len0 < 1e-6 || len1 < 1e-6)
      continue;

    double cosTheta = vecDot(v0, v1) / (len0 * len1);
    cosTheta = (std::max)(-1.0, (std::min)(1.0, cosTheta));
    const double angleRad = std::acos(cosTheta);

    if (angleRad < thresholdRad) {
      outFeatureFlags[i] = 1;
    }
  }
}

} // anonymous namespace

//--------------------------
//---- PUBLIC METHODS
//--------------------------

// Directly rebuilds ioPath.graph from its existing vertex order. The
// spur-pruning/gap-welding/chain-selection this used to do now happens once
// upstream, during extraction (zSlicer::_extractTrackPaths), which also
// already splits inner/outer tracks into separate PrintPaths - so by the
// time a path gets here it is already a single clean, ordered chain/loop;
// this step just validates it and rebuilds the graph object from scratch.
ZSPACE_INLINE bool zPrintSyn::buildSequencePrintPath(PrintPath &ioPath,
                                                     double weldTolerance,
                                                     double spurLengthThreshold,
                                                     std::string *outReason) {
  (void)weldTolerance;
  (void)spurLengthThreshold;

  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  zIntArray edgeConnects;
  fnGraph.getVertexPositions(positions);
  fnGraph.getEdgeData(edgeConnects);

  const int numV = static_cast<int>(positions.size());
  const int numE = static_cast<int>(edgeConnects.size() / 2);
  if (numV < 3 || numE < 1) {
    if (outReason) {
      std::ostringstream oss;
      oss << "Insufficient elements (vertices=" << numV << ", edges=" << numE
          << ")";
      *outReason = oss.str();
    }
    return false;
  }

  fnGraph.clear();
  fnGraph.create(positions, edgeConnects);
  return true;
}

// Fixed-spacing arc-length resample of ioPath.graph, carrying
// widths/heights/normals along. The sampling phase resets at each
// preserved feature/corner vertex (and at the open path's true endpoint)
// so those land exactly on an output point instead of being smoothed
// over; everywhere else, output points are placed at a uniform
// sampleLength apart regardless of the input's original vertex spacing
// (unlike a simple "subdivide long segments" pass, this can also make
// an already-dense input coarser).
ZSPACE_INLINE bool zPrintSyn::resamplePrintPath(PrintPath &ioPath) {
  zFnGraph fnGraph(ioPath.graph);
  const int N = fnGraph.numVertices();
  const float sampleLength =
      (i_Settings.ResampleLength > 1e-4f) ? i_Settings.ResampleLength : 0.005f;
  if (N < 2 || sampleLength <= 1e-6f)
    return false;

  const bool closed = isClosedGraph(ioPath.graph);

  zPointArray positions;
  fnGraph.getVertexPositions(positions);

  std::vector<int> sequence(N);
  for (int i = 0; i < N; i++)
    sequence[i] = i;

  std::vector<int> featureFlags;
  _detectFeatureAngles(positions, sequence, i_Settings.FeatureAngleTol,
                      featureFlags);

  const zFloatArray &weights = ioPath.widths;
  const zFloatArray &heights = ioPath.heights;
  const zVectorArray &normals = ioPath.normals;
  const bool hasHeights = !heights.empty();
  const bool hasNormals = !normals.empty();

  zPointArray outPoints;
  zFloatArray outWeights;
  zFloatArray outHeights;
  zVectorArray outNormals;
  zIntArray outEdges;

  auto appendPoint = [&](const zPoint &p, float w, float h,
                         const zVector &nrm) {
    if (!outPoints.empty()) {
      zVector d = vecSub(p, outPoints.back());
      if (d.length() < 1e-5)
        return;
    }
    outPoints.push_back(p);
    outWeights.push_back(w);
    if (hasHeights)
      outHeights.push_back(h);
    if (hasNormals)
      outNormals.push_back(nrm);
  };

  auto weightAt = [&](int id) {
    return (id < static_cast<int>(weights.size())) ? weights[id] : 0.04f;
  };
  auto heightAt = [&](int id) {
    return (hasHeights && id < static_cast<int>(heights.size()))
               ? heights[id]
               : 0.0f;
  };
  auto normalAt = [&](int id) {
    return (hasNormals && id < static_cast<int>(normals.size()))
               ? normals[id]
               : zVector(0, 0, 1);
  };

  const int numSegments = closed ? N : (N - 1);

  // Emit the path's start vertex and reset the sampling phase there.
  const int firstId = sequence[0];
  appendPoint(positions[firstId], weightAt(firstId), heightAt(firstId),
             normalAt(firstId));
  double distSinceSample = 0.0;

  for (int seg = 0; seg < numSegments; seg++) {
    int i0 = sequence[seg];
    int i1 = sequence[(seg + 1) % N];
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    float w0 = weightAt(i0), w1 = weightAt(i1);
    float h0 = heightAt(i0), h1 = heightAt(i1);
    zVector n0 = normalAt(i0), n1 = normalAt(i1);

    zVector segVec = vecSub(p1, p0);
    const double segLen = segVec.length();

    if (segLen > 1e-9) {
      // March along this segment, continuing the running phase carried
      // over from previous segments, placing one sample every
      // sampleLength regardless of how long/short this segment is.
      double remaining = sampleLength - distSinceSample;
      while (remaining < segLen) {
        const double t = remaining / segLen;
        zPoint interpP(p0.x + t * segVec.x, p0.y + t * segVec.y,
                       p0.z + t * segVec.z);
        float interpW = static_cast<float>(w0 + t * (w1 - w0));
        float interpH = static_cast<float>(h0 + t * (h1 - h0));
        zVector interpN(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
                        n0.z + t * (n1.z - n0.z));
        if (interpN.length() > 1e-6)
          interpN.normalize();
        else
          interpN = zVector(0, 0, 1);
        appendPoint(interpP, interpW, interpH, interpN);
        remaining += sampleLength;
      }
      distSinceSample = segLen - (remaining - sampleLength);
    }

    // Force-emit i1 (and reset the phase there) if it's the open path's
    // true endpoint or a preserved feature/corner vertex.
    const bool isLastOpenVertex = (!closed && seg == numSegments - 1);
    const bool isFeature =
        (i1 < int(featureFlags.size())) && featureFlags[i1] != 0;
    if (isLastOpenVertex || isFeature) {
      appendPoint(p1, w1, h1, n1);
      distSinceSample = 0.0;
    }
  }

  const int numPts = static_cast<int>(outPoints.size());
  if (numPts >= 2) {
    for (int i = 0; i + 1 < numPts; i++) {
      outEdges.push_back(i);
      outEdges.push_back(i + 1);
    }
    if (closed && numPts >= 3) {
      zVector d = vecSub(outPoints.back(), outPoints.front());
      if (d.length() < 1e-4) {
        outPoints.pop_back();
        outWeights.pop_back();
        if (hasHeights)
          outHeights.pop_back();
        if (hasNormals)
          outNormals.pop_back();
        outEdges.pop_back();
        outEdges.pop_back();
      }
      const int finalN = static_cast<int>(outPoints.size());
      if (finalN >= 3) {
        outEdges.push_back(finalN - 1);
        outEdges.push_back(0);
      }
    }
  }

  if (outPoints.size() < 3 || outEdges.size() < 4)
    return false;

  fnGraph.clear();
  fnGraph.create(outPoints, outEdges);
  ioPath.widths = outWeights;
  if (hasHeights)
    ioPath.heights = outHeights;
  if (hasNormals)
    ioPath.normals = outNormals;
  return true;
}

// Smooths localized width glitches (isolated spikes/dips) in ioPath.widths.
ZSPACE_INLINE void zPrintSyn::filterGlitches(PrintPath &ioPath) {
  const bool closed = isClosedGraph(ioPath.graph);
  zFloatArray &widths = ioPath.widths;
  const int n = static_cast<int>(widths.size());
  if (n <= 2)
    return;
  for (int i = 0; i < n; i++) {
    int prevId = (i == 0) ? (closed ? n - 1 : 0) : i - 1;
    int nextId = (i == n - 1) ? (closed ? 0 : n - 1) : i + 1;
    if (prevId == i || nextId == i)
      continue;
    float prevWidth = widths[prevId];
    float currentWidth = widths[i];
    float nextWidth = widths[nextId];
    if (prevWidth <= 0.0f || nextWidth <= 0.0f)
      continue;
    if (std::fabs(prevWidth - nextWidth) > 0.001f)
      continue;
    if (std::fabs(currentWidth - prevWidth) <= 0.001f)
      continue;
    widths[i] = prevWidth;
  }
}

} // namespace zSpace
