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
#include <limits>
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

// Signed area (shoelace formula) of a closed 2D polyline in the XY plane;
// positive = counter-clockwise, negative = clockwise.
double _signedArea2D(const zPointArray &positions) {
  double area = 0.0;
  const int n = static_cast<int>(positions.size());
  for (int i = 0; i < n; i++) {
    const zPoint &p0 = positions[i];
    const zPoint &p1 = positions[(i + 1) % n];
    area += (p0.x * p1.y - p1.x * p0.y);
  }
  return 0.5 * area;
}

// Reverses ioPath's point order in place (graph + widths/heights/normals),
// rebuilding the graph with the same open/closed edge structure it had.
void _reversePrintPath(PrintPath &ioPath) {
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int numV = static_cast<int>(positions.size());
  if (numV < 2)
    return;

  const bool closed = isClosedGraph(ioPath.graph);

  std::reverse(positions.begin(), positions.end());
  if (!ioPath.widths.empty())
    std::reverse(ioPath.widths.begin(), ioPath.widths.end());
  if (!ioPath.heights.empty())
    std::reverse(ioPath.heights.begin(), ioPath.heights.end());
  if (!ioPath.normals.empty())
    std::reverse(ioPath.normals.begin(), ioPath.normals.end());
  if (!ioPath.trackIds.empty())
    std::reverse(ioPath.trackIds.begin(), ioPath.trackIds.end());

  zIntArray newEdges;
  newEdges.reserve(numV * 2);
  for (int i = 0; i + 1 < numV; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }
  if (closed && numV >= 3) {
    newEdges.push_back(numV - 1);
    newEdges.push_back(0);
  }

  fnGraph.clear();
  fnGraph.create(positions, newEdges);
}

// Cumulative arc length walking a closed loop's positions (size N+1,
// including the wraparound segment back to index 0).
std::vector<double> _cumulativeArcLength(const zPointArray &positions) {
  const int n = static_cast<int>(positions.size());
  std::vector<double> cumLen(n + 1, 0.0);
  for (int i = 0; i < n; i++) {
    cumLen[i + 1] =
        cumLen[i] + vecSub(positions[(i + 1) % n], positions[i]).length();
  }
  return cumLen;
}

// Inserts an interpolated vertex at arc-length position s into ioPath's
// closed loop (splitting whichever segment contains s, snapping to an
// existing vertex instead if s lands within weldEps of one), then rotates
// the loop so that vertex becomes index 0. widths/heights/normals/
// trackIds are interpolated/rotated alongside (trackIds thresholded to
// 0/1, same convention as extraction).
bool _insertAndRotateToStart(PrintPath &ioPath, double s,
                             double weldEps = 1e-6) {
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 3 || !isClosedGraph(ioPath.graph))
    return false;

  const std::vector<double> cumLen = _cumulativeArcLength(positions);
  const double L = cumLen[n];
  if (L <= 1e-9)
    return false;
  s = std::fmod(s, L);
  if (s < 0.0)
    s += L;

  int seg = 0;
  for (; seg < n; seg++) {
    if (s <= cumLen[seg + 1] || seg == n - 1)
      break;
  }
  const int i0 = seg, i1 = (seg + 1) % n;
  const double segStart = cumLen[seg], segLen = cumLen[seg + 1] - segStart;
  const double t = (segLen > 1e-9)
                       ? (std::max)(0.0, (std::min)(1.0,
                                                    (s - segStart) / segLen))
                       : 0.0;

  const bool hasHeights = !ioPath.heights.empty();
  const bool hasNormals = !ioPath.normals.empty();
  const bool hasTrackIds = !ioPath.trackIds.empty();

  zPointArray newPositions = positions;
  zFloatArray newWidths = ioPath.widths;
  zFloatArray newHeights = ioPath.heights;
  zVectorArray newNormals = ioPath.normals;
  zFloatArray newTrackIds = ioPath.trackIds;
  int startIdx;

  const double tEps = weldEps / (std::max)(segLen, 1e-9);
  if (t <= tEps) {
    startIdx = i0;
  } else if (t >= 1.0 - tEps) {
    startIdx = i1;
  } else {
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    zPoint newP(p0.x + t * (p1.x - p0.x), p0.y + t * (p1.y - p0.y),
               p0.z + t * (p1.z - p0.z));

    const float w0 =
        (i0 < int(ioPath.widths.size())) ? ioPath.widths[i0] : 0.04f;
    const float w1 =
        (i1 < int(ioPath.widths.size())) ? ioPath.widths[i1] : 0.04f;
    const float newW = static_cast<float>(w0 + t * (w1 - w0));

    float newH = 0.0f;
    if (hasHeights) {
      const float h0 =
          (i0 < int(ioPath.heights.size())) ? ioPath.heights[i0] : 0.0f;
      const float h1 =
          (i1 < int(ioPath.heights.size())) ? ioPath.heights[i1] : 0.0f;
      newH = static_cast<float>(h0 + t * (h1 - h0));
    }
    zVector newN(0, 0, 1);
    if (hasNormals) {
      const zVector n0 = (i0 < int(ioPath.normals.size()))
                             ? ioPath.normals[i0]
                             : zVector(0, 0, 1);
      const zVector n1 = (i1 < int(ioPath.normals.size()))
                             ? ioPath.normals[i1]
                             : zVector(0, 0, 1);
      zVector n(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
               n0.z + t * (n1.z - n0.z));
      if (n.length() > 1e-6)
        n.normalize();
      else
        n = zVector(0, 0, 1);
      newN = n;
    }
    float newTrk = 0.0f;
    if (hasTrackIds) {
      const float trk0 =
          (i0 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i0] : 0.0f;
      const float trk1 =
          (i1 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i1] : 0.0f;
      newTrk = (trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f : 0.0f;
    }

    newPositions.insert(newPositions.begin() + i1, newP);
    newWidths.insert(newWidths.begin() + i1, newW);
    if (hasHeights)
      newHeights.insert(newHeights.begin() + i1, newH);
    if (hasNormals)
      newNormals.insert(newNormals.begin() + i1, newN);
    if (hasTrackIds)
      newTrackIds.insert(newTrackIds.begin() + i1, newTrk);
    startIdx = i1;
  }

  std::rotate(newPositions.begin(), newPositions.begin() + startIdx,
             newPositions.end());
  std::rotate(newWidths.begin(), newWidths.begin() + startIdx,
             newWidths.end());
  if (hasHeights)
    std::rotate(newHeights.begin(), newHeights.begin() + startIdx,
               newHeights.end());
  if (hasNormals)
    std::rotate(newNormals.begin(), newNormals.begin() + startIdx,
               newNormals.end());
  if (hasTrackIds)
    std::rotate(newTrackIds.begin(), newTrackIds.begin() + startIdx,
               newTrackIds.end());

  const int newN = static_cast<int>(newPositions.size());
  zIntArray newEdges;
  newEdges.reserve(newN * 2);
  for (int i = 0; i + 1 < newN; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }
  newEdges.push_back(newN - 1);
  newEdges.push_back(0);

  fnGraph.clear();
  fnGraph.create(newPositions, newEdges);
  ioPath.widths = newWidths;
  if (hasHeights)
    ioPath.heights = newHeights;
  if (hasNormals)
    ioPath.normals = newNormals;
  if (hasTrackIds)
    ioPath.trackIds = newTrackIds;
  return true;
}

// Finds the arc-length position of the point on path's closed loop nearest
// to target, searching only segments whose both endpoints are inner-tagged
// (or every segment, if the path carries no track distinction at all) -
// the hard "stays on the inner ring" constraint. Returns false if no
// eligible segment exists.
bool _nearestInnerArcPos(const PrintPath &path, const zPoint &target,
                         double &outArcPos) {
  zFnGraph fnGraph(const_cast<zObjGraph &>(path.graph));
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 3 || !isClosedGraph(path.graph))
    return false;

  const bool hasTrackIds = !path.trackIds.empty();
  auto isOuter = [&](int id) {
    return hasTrackIds && id < int(path.trackIds.size()) &&
           path.trackIds[id] >= 0.5f;
  };

  const std::vector<double> cumLen = _cumulativeArcLength(positions);

  double bestDistSq = (std::numeric_limits<double>::max)();
  double bestArcPos = 0.0;
  bool found = false;

  for (int i = 0; i < n; i++) {
    const int i0 = i, i1 = (i + 1) % n;
    if (isOuter(i0) || isOuter(i1))
      continue; // keep the search strictly inside the inner segment

    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    zVector seg = vecSub(p1, p0);
    const double segLenSq = vecDot(seg, seg);
    double t = 0.0;
    if (segLenSq > 1e-12) {
      t = vecDot(vecSub(target, p0), seg) / segLenSq;
      t = (std::max)(0.0, (std::min)(1.0, t));
    }
    zPoint proj(p0.x + t * seg.x, p0.y + t * seg.y, p0.z + t * seg.z);
    const double distSq = vecDistSq(target, proj);
    if (distSq < bestDistSq) {
      bestDistSq = distSq;
      bestArcPos = cumLen[i] + t * (cumLen[i + 1] - cumLen[i]);
      found = true;
    }
  }

  if (found)
    outArcPos = bestArcPos;
  return found;
}

// Drops a closed loop's wraparound edge (last vertex -> first vertex), so
// the seam sits exactly at the loop's current start/end instead of being
// welded shut. Point order and count are unchanged - only the closing
// edge goes away. No-op if the path isn't closed.
void _openLoopAtStart(PrintPath &ioPath) {
  if (!isClosedGraph(ioPath.graph))
    return;
  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 2)
    return;

  zIntArray newEdges;
  newEdges.reserve((n - 1) * 2);
  for (int i = 0; i + 1 < n; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }

  fnGraph.clear();
  fnGraph.create(positions, newEdges);
}

// Arc length of positions' inner-tagged prefix (index 0 is always inner
// and inner is a contiguous prefix, per _extractTrackPaths), or the whole
// loop's length if trackIds carries no track distinction at all.
double _innerSegmentLength(const zPointArray &positions,
                           const zFloatArray &trackIds,
                           const std::vector<double> &cumLen) {
  const int n = static_cast<int>(positions.size());
  int innerEnd = n;
  if (!trackIds.empty()) {
    innerEnd = 0;
    while (innerEnd < n && innerEnd < int(trackIds.size()) &&
           trackIds[innerEnd] < 0.5f)
      innerEnd++;
  }
  return cumLen[innerEnd];
}

// Cuts the closed ring open right before its (already finalized) start
// point: walks backward seamSize along the loop from the start to find a
// new end point (inserting an interpolated vertex there, unless it snaps
// onto an existing one), keeps start..that new end, and drops the small
// remaining arc between the new end and the start. seamSize <= 0
// degenerates to _openLoopAtStart's plain "drop the closing edge" (no new
// point - the seed of this whole function is a seamSize of exactly 0).
void _cutGapBeforeStart(PrintPath &ioPath, double seamSize) {
  if (seamSize <= 1e-9) {
    _openLoopAtStart(ioPath);
    return;
  }
  if (!isClosedGraph(ioPath.graph))
    return;

  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int n = static_cast<int>(positions.size());
  if (n < 3)
    return;

  const std::vector<double> cumLen = _cumulativeArcLength(positions);
  const double L = cumLen[n];
  if (L <= 1e-9)
    return;

  const double gap = (std::min)(seamSize, L * 0.5); // never eat more than half the loop
  const double targetArcPos = std::fmod(L - gap, L);

  int seg = n - 1;
  for (int i = 0; i < n; i++) {
    if (targetArcPos <= cumLen[i + 1]) {
      seg = i;
      break;
    }
  }
  const int i0 = seg, i1 = (seg + 1) % n;
  const double segStart = cumLen[seg], segLen = cumLen[seg + 1] - segStart;
  const double t = (segLen > 1e-9)
                       ? (std::max)(0.0, (std::min)(1.0,
                                                    (targetArcPos - segStart) /
                                                        segLen))
                       : 0.0;

  const bool hasHeights = !ioPath.heights.empty();
  const bool hasNormals = !ioPath.normals.empty();
  const bool hasTrackIds = !ioPath.trackIds.empty();

  zPointArray keptPositions;
  zFloatArray keptWidths;
  zFloatArray keptHeights;
  zVectorArray keptNormals;
  zFloatArray keptTrackIds;

  auto appendIdx = [&](int id) {
    keptPositions.push_back(positions[id]);
    keptWidths.push_back(id < int(ioPath.widths.size()) ? ioPath.widths[id]
                                                        : 0.04f);
    if (hasHeights)
      keptHeights.push_back(id < int(ioPath.heights.size())
                                ? ioPath.heights[id]
                                : 0.0f);
    if (hasNormals)
      keptNormals.push_back(id < int(ioPath.normals.size())
                                ? ioPath.normals[id]
                                : zVector(0, 0, 1));
    if (hasTrackIds)
      keptTrackIds.push_back(id < int(ioPath.trackIds.size())
                                 ? ioPath.trackIds[id]
                                 : 0.0f);
  };

  // Keep the surviving arc from the start (index 0) through i0, then
  // append the new end point (interpolated, or i1/i0 if it snaps onto an
  // existing vertex) - producing an open chain from start to this new
  // end, with the small [end, start) gap dropped.
  for (int id = 0; id <= i0; id++)
    appendIdx(id);

  const double tEps = 1e-6 / (std::max)(segLen, 1e-9);
  if (t >= 1.0 - tEps) {
    if (i1 != i0)
      appendIdx(i1);
  } else if (t > tEps) {
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    zPoint newP(p0.x + t * (p1.x - p0.x), p0.y + t * (p1.y - p0.y),
               p0.z + t * (p1.z - p0.z));
    const float w0 =
        (i0 < int(ioPath.widths.size())) ? ioPath.widths[i0] : 0.04f;
    const float w1 =
        (i1 < int(ioPath.widths.size())) ? ioPath.widths[i1] : 0.04f;
    keptPositions.push_back(newP);
    keptWidths.push_back(static_cast<float>(w0 + t * (w1 - w0)));
    if (hasHeights) {
      const float h0 =
          (i0 < int(ioPath.heights.size())) ? ioPath.heights[i0] : 0.0f;
      const float h1 =
          (i1 < int(ioPath.heights.size())) ? ioPath.heights[i1] : 0.0f;
      keptHeights.push_back(static_cast<float>(h0 + t * (h1 - h0)));
    }
    if (hasNormals) {
      const zVector n0 = (i0 < int(ioPath.normals.size()))
                             ? ioPath.normals[i0]
                             : zVector(0, 0, 1);
      const zVector n1 = (i1 < int(ioPath.normals.size()))
                             ? ioPath.normals[i1]
                             : zVector(0, 0, 1);
      zVector nrm(n0.x + t * (n1.x - n0.x), n0.y + t * (n1.y - n0.y),
                 n0.z + t * (n1.z - n0.z));
      if (nrm.length() > 1e-6)
        nrm.normalize();
      else
        nrm = zVector(0, 0, 1);
      keptNormals.push_back(nrm);
    }
    if (hasTrackIds) {
      const float trk0 =
          (i0 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i0] : 0.0f;
      const float trk1 =
          (i1 < int(ioPath.trackIds.size())) ? ioPath.trackIds[i1] : 0.0f;
      keptTrackIds.push_back((trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f
                                                                 : 0.0f);
    }
  }
  // t <= tEps: targetArcPos already snapped onto i0, which is already the
  // last appended point - nothing more to add.

  const int keptN = static_cast<int>(keptPositions.size());
  if (keptN < 2)
    return;

  zIntArray newEdges;
  newEdges.reserve((keptN - 1) * 2);
  for (int i = 0; i + 1 < keptN; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }

  fnGraph.clear();
  fnGraph.create(keptPositions, newEdges);
  ioPath.widths = keptWidths;
  if (hasHeights)
    ioPath.heights = keptHeights;
  if (hasNormals)
    ioPath.normals = keptNormals;
  if (hasTrackIds)
    ioPath.trackIds = keptTrackIds;
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
  const zFloatArray &trackIds = ioPath.trackIds;
  const bool hasHeights = !heights.empty();
  const bool hasNormals = !normals.empty();
  const bool hasTrackIds = !trackIds.empty();

  zPointArray outPoints;
  zFloatArray outWeights;
  zFloatArray outHeights;
  zVectorArray outNormals;
  zFloatArray outTrackIds;
  zIntArray outEdges;

  auto appendPoint = [&](const zPoint &p, float w, float h,
                         const zVector &nrm, float trk) {
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
    if (hasTrackIds)
      outTrackIds.push_back(trk);
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
  auto trackAt = [&](int id) {
    return (hasTrackIds && id < static_cast<int>(trackIds.size()))
               ? trackIds[id]
               : 0.0f;
  };

  const int numSegments = closed ? N : (N - 1);

  // Emit the path's start vertex and reset the sampling phase there.
  const int firstId = sequence[0];
  appendPoint(positions[firstId], weightAt(firstId), heightAt(firstId),
             normalAt(firstId), trackAt(firstId));
  double distSinceSample = 0.0;

  for (int seg = 0; seg < numSegments; seg++) {
    int i0 = sequence[seg];
    int i1 = sequence[(seg + 1) % N];
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    float w0 = weightAt(i0), w1 = weightAt(i1);
    float h0 = heightAt(i0), h1 = heightAt(i1);
    zVector n0 = normalAt(i0), n1 = normalAt(i1);
    float trk0 = trackAt(i0), trk1 = trackAt(i1);

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
        const float interpTrk = (trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f
                                                                    : 0.0f;
        appendPoint(interpP, interpW, interpH, interpN, interpTrk);
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
      appendPoint(p1, w1, h1, n1, trk1);
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
        if (hasTrackIds)
          outTrackIds.pop_back();
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
  if (hasTrackIds)
    ioPath.trackIds = outTrackIds;
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

// Orients one layer's ring path CCW in place (extraction already produces
// at most one ring per layer, already inner-before-outer ordered).
ZSPACE_INLINE bool
zPrintSyn::orientAndWeldLayer(std::vector<PrintPath> &ringPaths) {
  if (ringPaths.empty())
    return false;

  PrintPath &path = ringPaths[0];
  zFnGraph fnGraph(path.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  if (positions.size() < 3)
    return false;
  if (_signedArea2D(positions) < 0.0) {
    _reversePrintPath(path);
  }
  return true;
}

// Aligns every layer's ring start point to the same barycentric point on
// the underlying mesh, constrained to the inner track and nudged by
// seamOffset. Seeds from the first usable layer by walking walkDistance
// from its extracted start along the inner segment (clamped so it can't
// cross into outer), parameterizes that point against the seed layer's 2D
// section mesh, then realizes the same parametric point on every other
// layer's mesh and rotates that layer's ring to start at its nearest
// inner-segment point - tracking the same physical spot on the model
// rather than a fixed 2D region, so it stays close layer-to-layer even
// through taper/twist. Finally cuts every ring open at its (now aligned)
// start via _cutGapBeforeStart(seamSize), so each layer ends as an open
// chain ready for chainLayersIntoOnePath.
ZSPACE_INLINE bool zPrintSyn::alignLayerStartPoints(
    std::vector<std::vector<PrintPath>> &ioPaths,
    const zObjMeshArray &sectionMeshes2D, double walkDistance,
    double seamOffset, double seamSize) {
  const int totalLayers = int(ioPaths.size());
  if (totalLayers == 0 || int(sectionMeshes2D.size()) < totalLayers)
    return false;

  int seedLayer = -1;
  for (int l = 0; l < totalLayers; l++) {
    if (!ioPaths[l].empty()) {
      zFnGraph fnGraph(ioPaths[l][0].graph);
      if (fnGraph.numVertices() >= 3) {
        seedLayer = l;
        break;
      }
    }
  }
  if (seedLayer < 0)
    return false;

  // 1. Seed: walk walkDistance + seamOffset from the current start along
  // the inner segment (index 0 is always inner, as a contiguous prefix,
  // per _extractTrackPaths), clamped so it can't cross into outer.
  PrintPath &seedPath = ioPaths[seedLayer][0];
  {
    zFnGraph fnGraph(seedPath.graph);
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    const std::vector<double> cumLen = _cumulativeArcLength(positions);
    const double innerLen =
        _innerSegmentLength(positions, seedPath.trackIds, cumLen);
    const double clampedWalk = (std::max)(
        0.0,
        (std::min)(walkDistance + seamOffset, (std::max)(0.0, innerLen - 1e-4)));
    if (!_insertAndRotateToStart(seedPath, clampedWalk))
      return false;
  }

  zPointArray seedPts;
  {
    zFnGraph fnGraph(seedPath.graph);
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    if (positions.empty())
      return false;
    seedPts.push_back(positions[0]);
    seedPts.push_back(positions[0]); // duplicate: parameterizeGraph needs an edge
  }

  zObjGraph seedGraph;
  zIntArray seedEdges = {0, 1};
  zFnGraph(seedGraph).create(seedPts, seedEdges);

  ParameterGraph seedParams;
  if (!parameterizeGraph(seedGraph,
                         const_cast<zObjMesh &>(sectionMeshes2D[seedLayer]),
                         seedParams))
    return false;

  // 2. Every layer: realize the same mesh-parametric point on its own 2D
  // section mesh, then rotate its ring to start at the nearest inner-
  // segment point to it. No sequential dependency between layers (the
  // target is the fixed seedParams, not the previous layer's result), so
  // this is safe to run in parallel.
  int alignedCount = 0;
#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    if (ioPaths[l].empty())
      continue;
    PrintPath &path = ioPaths[l][0];
    zFnGraph fnGraph(path.graph);
    if (fnGraph.numVertices() < 3)
      continue;

    zObjGraph realized;
    if (!realizeParameterGraph(seedParams,
                               const_cast<zObjMesh &>(sectionMeshes2D[l]),
                               realized))
      continue;
    zFnGraph fnRealized(realized);
    zPointArray realizedPts;
    fnRealized.getVertexPositions(realizedPts);
    if (realizedPts.empty())
      continue;

    double arcPos = 0.0;
    if (!_nearestInnerArcPos(path, realizedPts[0], arcPos))
      continue;

    zPointArray positionsL;
    fnGraph.getVertexPositions(positionsL);
    const std::vector<double> cumLenL = _cumulativeArcLength(positionsL);
    const double innerLenL =
        _innerSegmentLength(positionsL, path.trackIds, cumLenL);
    arcPos = (std::max)(
        0.0, (std::min)(arcPos + seamOffset, (std::max)(0.0, innerLenL - 1e-4)));

    if (_insertAndRotateToStart(path, arcPos)) {
#pragma omp atomic
      alignedCount++;
    }
  }

  // 3. Cut every layer's ring open right before its (now aligned) start
  // point, leaving a seamSize gap between the new end point and the
  // start - ready to be chained layer-to-layer.
  for (int l = 0; l < totalLayers; l++) {
    if (ioPaths[l].empty())
      continue;
    _cutGapBeforeStart(ioPaths[l][0], seamSize);
  }

  return alignedCount > 0;
}

// Concatenates every layer's single (post-alignLayerStartPoints, open)
// PrintPath end-to-start into one continuous bottom-to-top PrintPath.
// Simple sequential concatenation already connects each layer's last
// point to the next layer's first point via the edge between them.
ZSPACE_INLINE bool zPrintSyn::chainLayersIntoOnePath(
    const std::vector<std::vector<PrintPath>> &layeredPaths,
    PrintPath &outCombined) {
  outCombined = PrintPath();

  bool hasHeights = false, hasNormals = false, hasTrackIds = false;
  for (const auto &ringPaths : layeredPaths) {
    if (!ringPaths.empty()) {
      if (!ringPaths[0].heights.empty())
        hasHeights = true;
      if (!ringPaths[0].normals.empty())
        hasNormals = true;
      if (!ringPaths[0].trackIds.empty())
        hasTrackIds = true;
    }
  }

  zPointArray combinedPositions;
  zFloatArray combinedWidths;
  zFloatArray combinedHeights;
  zVectorArray combinedNormals;
  zFloatArray combinedTrackIds;

  for (const auto &ringPaths : layeredPaths) {
    if (ringPaths.empty())
      continue;
    const PrintPath &layerPath = ringPaths[0]; // one path per layer, post
                                               // alignLayerStartPoints

    zFnGraph fnGraph(const_cast<zObjGraph &>(layerPath.graph));
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    if (positions.empty())
      continue;

    combinedPositions.insert(combinedPositions.end(), positions.begin(),
                             positions.end());
    combinedWidths.insert(combinedWidths.end(), layerPath.widths.begin(),
                          layerPath.widths.end());
    if (hasHeights) {
      if (!layerPath.heights.empty())
        combinedHeights.insert(combinedHeights.end(),
                               layerPath.heights.begin(),
                               layerPath.heights.end());
      else
        combinedHeights.insert(combinedHeights.end(), positions.size(), 0.0f);
    }
    if (hasNormals) {
      if (!layerPath.normals.empty())
        combinedNormals.insert(combinedNormals.end(),
                               layerPath.normals.begin(),
                               layerPath.normals.end());
      else
        combinedNormals.insert(combinedNormals.end(), positions.size(),
                               zVector(0, 0, 1));
    }
    if (hasTrackIds) {
      if (!layerPath.trackIds.empty())
        combinedTrackIds.insert(combinedTrackIds.end(),
                                layerPath.trackIds.begin(),
                                layerPath.trackIds.end());
      else
        combinedTrackIds.insert(combinedTrackIds.end(), positions.size(),
                                0.0f);
    }
  }

  const int total = static_cast<int>(combinedPositions.size());
  if (total < 2)
    return false;

  zIntArray combinedEdges;
  combinedEdges.reserve((total - 1) * 2);
  for (int i = 0; i + 1 < total; i++) {
    combinedEdges.push_back(i);
    combinedEdges.push_back(i + 1);
  }

  zFnGraph fnOut(outCombined.graph);
  fnOut.create(combinedPositions, combinedEdges);
  outCombined.widths = combinedWidths;
  if (hasHeights)
    outCombined.heights = combinedHeights;
  if (hasNormals)
    outCombined.normals = combinedNormals;
  if (hasTrackIds)
    outCombined.trackIds = combinedTrackIds;
  return true;
}

} // namespace zSpace
