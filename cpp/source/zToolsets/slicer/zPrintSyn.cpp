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
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

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

    double angleRad;
    if (turnAngleAt(positions[prevId], positions[curId], positions[nextId],
                    angleRad) &&
        angleRad < thresholdRad) {
      outFeatureFlags[i] = 1;
    }
  }
}

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

std::vector<double> _cumulativeArcLength(const zPointArray &positions) {
  const int n = static_cast<int>(positions.size());
  std::vector<double> cumLen(n + 1, 0.0);
  for (int i = 0; i < n; i++) {
    cumLen[i + 1] =
        cumLen[i] + vecSub(positions[(i + 1) % n], positions[i]).length();
  }
  return cumLen;
}

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
  const double t =
      (segLen > 1e-9)
          ? (std::max)(0.0, (std::min)(1.0, (s - segStart) / segLen))
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
      const zVector n0 = (i0 < int(ioPath.normals.size())) ? ioPath.normals[i0]
                                                           : zVector(0, 0, 1);
      const zVector n1 = (i1 < int(ioPath.normals.size())) ? ioPath.normals[i1]
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
  std::rotate(newWidths.begin(), newWidths.begin() + startIdx, newWidths.end());
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

  const double gap =
      (std::min)(seamSize, L * 0.5); // never eat more than half the loop
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
  const double t =
      (segLen > 1e-9)
          ? (std::max)(0.0, (std::min)(1.0, (targetArcPos - segStart) / segLen))
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
      keptHeights.push_back(id < int(ioPath.heights.size()) ? ioPath.heights[id]
                                                            : 0.0f);
    if (hasNormals)
      keptNormals.push_back(id < int(ioPath.normals.size()) ? ioPath.normals[id]
                                                            : zVector(0, 0, 1));
    if (hasTrackIds)
      keptTrackIds.push_back(
          id < int(ioPath.trackIds.size()) ? ioPath.trackIds[id] : 0.0f);
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
      const zVector n0 = (i0 < int(ioPath.normals.size())) ? ioPath.normals[i0]
                                                           : zVector(0, 0, 1);
      const zVector n1 = (i1 < int(ioPath.normals.size())) ? ioPath.normals[i1]
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
      keptTrackIds.push_back((trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f : 0.0f);
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

// Extracts the zero-isocontour of pathVals as raw position/edge arrays
// (marching triangles), interpolating widths and the outer/inner track tag
// onto each new crossing point. Left unsplit/uncleaned - _extractTrackPaths
// does the spur-pruning, gap-welding and inner/outer split downstream.
// Moved here from zSlicer_Field.cpp - isocontour extraction is a printsyn
// (toolpath synthesis) responsibility, not a field-computation one.
void _extractRawIsoContour(
    zObjMeshScalarField &pathFieldMesh, const zScalarArray &pathVals,
    const zScalarArray &widthVals, const zScalarArray &trackIdVals,
    zPointArray &outPositions, zIntArray &outEdgeConnects,
    zFloatArray &outVertexWeights, zFloatArray &outVertexTrackIds) {
  outPositions.clear();
  outEdgeConnects.clear();
  outVertexWeights.clear();
  outVertexTrackIds.clear();

  zFnMesh fnMesh(pathFieldMesh);
  const zPoint *rawPositions = fnMesh.getRawVertexPositions();
  const int numVerts = fnMesh.numVertices();
  const int numFaces = fnMesh.numPolygons();
  const bool hasTrackIds = trackIdVals.size() == static_cast<size_t>(numVerts);
  if (!rawPositions || numVerts == 0 || numFaces == 0 ||
      pathVals.size() != static_cast<size_t>(numVerts) ||
      widthVals.size() != static_cast<size_t>(numVerts)) {
    return;
  }

  zPointArray &isoPositions = outPositions;
  zIntArray &isoEdgeConnects = outEdgeConnects;

  std::unordered_map<unsigned long long, int> edgeCrossMap;
  edgeCrossMap.reserve(numFaces / 4);

  auto getEdgeIntersection = [&](int v0, int v1) -> int {
    const unsigned long long edgeKey = makeEdgeKey(v0, v1);

    auto it = edgeCrossMap.find(edgeKey);
    if (it != edgeCrossMap.end())
      return it->second;

    const double s0 = pathVals[v0];
    const double s1 = pathVals[v1];
    const double denom = s1 - s0;
    const double t = (std::abs(denom) > 1e-12)
                         ? (std::max)(0.0, (std::min)(1.0, -s0 / denom))
                         : 0.5;

    const zPoint &p0 = rawPositions[v0], &p1 = rawPositions[v1];
    isoPositions.push_back(zPoint(p0.x + t * (p1.x - p0.x),
                                  p0.y + t * (p1.y - p0.y),
                                  p0.z + t * (p1.z - p0.z)));
    outVertexWeights.push_back(static_cast<float>(
        widthVals[v0] + t * (widthVals[v1] - widthVals[v0])));

    if (hasTrackIds) {
      const double id0 = trackIdVals[v0];
      const double id1 = trackIdVals[v1];
      const double idInterp = id0 + t * (id1 - id0);
      outVertexTrackIds.push_back(idInterp >= 0.5 ? 1.0f : 0.0f);
    } else {
      outVertexTrackIds.push_back(0.0f);
    }

    const int newIdx = static_cast<int>(isoPositions.size()) - 1;
    edgeCrossMap[edgeKey] = newIdx;
    return newIdx;
  };

  zIntArray polyConnects, polyCounts;
  fnMesh.getPolygonData(polyConnects, polyCounts);

  size_t cOffset = 0;
  for (int count : polyCounts) {
    for (int tri = 1; tri + 1 < count; tri++) {
      const int v0 = polyConnects[cOffset];
      const int v1 = polyConnects[cOffset + tri];
      const int v2 = polyConnects[cOffset + tri + 1];

      const double s0 = pathVals[v0], s1 = pathVals[v1], s2 = pathVals[v2];
      const int mask = (s0 > 0.0) | ((s1 > 0.0) << 1) | ((s2 > 0.0) << 2);
      if (mask == 0 || mask == 7)
        continue;

      int p0 = -1, p1 = -1;
      switch (mask) {
      case 1:
      case 6:
        p0 = getEdgeIntersection(v0, v1);
        p1 = getEdgeIntersection(v2, v0);
        break;
      case 2:
      case 5:
        p0 = getEdgeIntersection(v0, v1);
        p1 = getEdgeIntersection(v1, v2);
        break;
      case 3:
      case 4:
        p0 = getEdgeIntersection(v1, v2);
        p1 = getEdgeIntersection(v2, v0);
        break;
      }

      if (p0 != p1 && p0 >= 0 && p1 >= 0) {
        isoEdgeConnects.push_back(p0);
        isoEdgeConnects.push_back(p1);
      }
    }
    cOffset += count;
  }
}

// Spur-pruning + gap-welding + chain/loop extraction: splits the raw
// iso-contour, walks it into a single oriented, inner-before-outer closed
// ring (plus bracing open chains), and merges every surviving loop/chain
// into one PrintPath.
void _extractTrackPaths(const zPointArray &rawPositions,
                        const zIntArray &rawEdges, const zFloatArray &rawWidths,
                        const zFloatArray &rawTrackIds, double weldTolerance,
                        double spurLengthThreshold,
                        std::vector<PrintPath> &outPaths) {
  outPaths.clear();
  const int numV = static_cast<int>(rawPositions.size());
  if (numV < 3 || rawEdges.size() < 6)
    return;

  const bool hasTrackIds = rawTrackIds.size() == static_cast<size_t>(numV);
  auto isOuterVertex = [&](int id) -> bool {
    return hasTrackIds && id < int(rawTrackIds.size()) &&
           rawTrackIds[id] >= 0.5f;
  };

  std::vector<std::vector<int>> adj(numV);
  std::unordered_set<unsigned long long> edgeSet;
  for (size_t e = 0; e + 1 < rawEdges.size(); e += 2) {
    int a = rawEdges[e], b = rawEdges[e + 1];
    if (a < 0 || b < 0 || a >= numV || b >= numV || a == b)
      continue;
    if (vecDistSq(rawPositions[b], rawPositions[a]) < 1e-12)
      continue;
    const unsigned long long key = makeEdgeKey(a, b);
    if (!edgeSet.insert(key).second)
      continue;
    adj[a].push_back(b);
    adj[b].push_back(a);
  }

  std::vector<std::unordered_set<int>> activeAdj(numV);
  for (int i = 0; i < numV; i++)
    for (int nb : adj[i])
      activeAdj[i].insert(nb);

  bool prunedAny = true;
  while (prunedAny) {
    prunedAny = false;
    for (int i = 0; i < numV; i++) {
      if (activeAdj[i].size() >= 3) {
        std::vector<int> neighbors(activeAdj[i].begin(), activeAdj[i].end());
        for (int startNb : neighbors) {
          std::vector<int> branchPath = {i, startNb};
          int prev = i;
          int curr = startNb;
          bool isDeadEnd = false;
          double branchLen =
              vecSub(rawPositions[curr], rawPositions[prev]).length();

          while (true) {
            if (activeAdj[curr].size() == 1) {
              isDeadEnd = true;
              break;
            }
            if (activeAdj[curr].size() != 2)
              break;
            int next = -1;
            for (int nb : activeAdj[curr]) {
              if (nb != prev) {
                next = nb;
                break;
              }
            }
            if (next < 0 || next == i)
              break;
            branchLen +=
                vecSub(rawPositions[next], rawPositions[curr]).length();
            prev = curr;
            curr = next;
            branchPath.push_back(curr);
          }

          if (isDeadEnd && branchLen < spurLengthThreshold) {
            for (size_t k = 0; k + 1 < branchPath.size(); k++) {
              int u = branchPath[k];
              int v = branchPath[k + 1];
              activeAdj[u].erase(v);
              activeAdj[v].erase(u);
            }
            prunedAny = true;
          }
        }
      }
    }
  }

  std::vector<std::vector<int>> openChains;
  std::vector<std::vector<int>> closedLoops;
  std::unordered_set<unsigned long long> visitedEdges;

  auto walkChain = [&](int start,
                       int next) -> std::pair<std::vector<int>, bool> {
    std::vector<int> seq = {start};
    int prev = start, curr = next;
    visitedEdges.insert(makeEdgeKey(prev, curr));

    while (true) {
      if (curr == start)
        return {seq, true};
      seq.push_back(curr);
      if (activeAdj[curr].size() != 2)
        break;

      int nxt = -1;
      for (int nb : activeAdj[curr]) {
        if (nb != prev) {
          nxt = nb;
          break;
        }
      }
      if (nxt < 0 || visitedEdges.count(makeEdgeKey(curr, nxt)))
        break;
      visitedEdges.insert(makeEdgeKey(curr, nxt));
      prev = curr;
      curr = nxt;
    }
    return {seq, false};
  };

  for (int i = 0; i < numV; i++) {
    if (activeAdj[i].size() == 1) {
      int next = *activeAdj[i].begin();
      if (visitedEdges.count(makeEdgeKey(i, next)))
        continue;
      openChains.push_back(walkChain(i, next).first);
    }
  }

  for (int i = 0; i < numV; i++) {
    if (activeAdj[i].size() == 2) {
      for (int nb : activeAdj[i]) {
        if (visitedEdges.count(makeEdgeKey(i, nb)))
          continue;

        auto [loop, closed] = walkChain(i, nb);
        if (closed && loop.size() >= 3) {
          closedLoops.push_back(std::move(loop));
        } else if (!loop.empty()) {
          openChains.push_back(std::move(loop));
        }
      }
    }
  }

  // 4. Self-close any open chain whose two ends are within weldTolerance -
  // safe unconditionally now that nothing is track-filtered, so the
  // boundary ring can never be mistakenly split by a filtering artifact
  // that this step would then need to leave alone.
  const double weldTolSq = weldTolerance * weldTolerance;
  for (auto it = openChains.begin(); it != openChains.end();) {
    if (it->size() >= 3 && vecDistSq(rawPositions[it->front()],
                                     rawPositions[it->back()]) <= weldTolSq) {
      closedLoops.push_back(*it);
      it = openChains.erase(it);
      continue;
    }
    ++it;
  }

  // 5. Weld remaining open-chain fragments (e.g. bracing ribs split at
  // triangle boundaries) end-to-end by nearest-endpoint matching.
  bool mergedAny = true;
  while (mergedAny && openChains.size() > 1) {
    mergedAny = false;
    size_t bestI = 0, bestJ = 0;
    int bestMode = -1;
    double bestDistSq = weldTolSq;

    for (size_t i = 0; i < openChains.size(); i++) {
      for (size_t j = i + 1; j < openChains.size(); j++) {
        const auto &c1 = openChains[i];
        const auto &c2 = openChains[j];

        double d0 =
            vecDistSq(rawPositions[c1.back()], rawPositions[c2.front()]);
        if (d0 < bestDistSq) {
          bestDistSq = d0;
          bestI = i;
          bestJ = j;
          bestMode = 0;
        }
        double d1 = vecDistSq(rawPositions[c1.back()], rawPositions[c2.back()]);
        if (d1 < bestDistSq) {
          bestDistSq = d1;
          bestI = i;
          bestJ = j;
          bestMode = 1;
        }
        double d2 =
            vecDistSq(rawPositions[c1.front()], rawPositions[c2.back()]);
        if (d2 < bestDistSq) {
          bestDistSq = d2;
          bestI = i;
          bestJ = j;
          bestMode = 2;
        }
        double d3 =
            vecDistSq(rawPositions[c1.front()], rawPositions[c2.front()]);
        if (d3 < bestDistSq) {
          bestDistSq = d3;
          bestI = i;
          bestJ = j;
          bestMode = 3;
        }
      }
    }

    if (bestMode >= 0) {
      std::vector<int> &c1 = openChains[bestI];
      std::vector<int> c2 = openChains[bestJ];

      if (bestMode == 0) {
        c1.insert(c1.end(), c2.begin(), c2.end());
      } else if (bestMode == 1) {
        std::reverse(c2.begin(), c2.end());
        c1.insert(c1.end(), c2.begin(), c2.end());
      } else if (bestMode == 2) {
        c2.insert(c2.end(), c1.begin(), c1.end());
        c1 = c2;
      } else if (bestMode == 3) {
        std::reverse(c2.begin(), c2.end());
        c2.insert(c2.end(), c1.begin(), c1.end());
        c1 = c2;
      }

      openChains.erase(openChains.begin() + bestJ);
      mergedAny = true;

      if (c1.size() >= 3 && vecDistSq(rawPositions[c1.front()],
                                      rawPositions[c1.back()]) <= weldTolSq) {
        closedLoops.push_back(c1);
        openChains.erase(openChains.begin() + bestI);
      }
    }
  }

  if (closedLoops.empty() && openChains.empty())
    return;

  // 6. Orient the longest closed loop (the boundary ring) CCW, and - if it
  // mixes inner/outer vertices - rotate it (it's already one continuous
  // loop, so no cutting/rewelding needed) so the walk starts at the
  // outer-to-inner transition and visits every inner vertex before any
  // outer one.
  if (!closedLoops.empty()) {
    size_t mainIdx = 0;
    for (size_t i = 1; i < closedLoops.size(); i++) {
      if (closedLoops[i].size() > closedLoops[mainIdx].size())
        mainIdx = i;
    }
    std::vector<int> &mainRing = closedLoops[mainIdx];
    const int m = static_cast<int>(mainRing.size());

    zPointArray ringPositions(m);
    for (int i = 0; i < m; i++)
      ringPositions[i] = rawPositions[mainRing[i]];
    if (signedArea2D(ringPositions) < 0.0)
      std::reverse(mainRing.begin(), mainRing.end());

    if (hasTrackIds) {
      bool hasInner = false, hasOuter = false;
      for (int id : mainRing) {
        if (isOuterVertex(id))
          hasOuter = true;
        else
          hasInner = true;
      }
      if (hasInner && hasOuter) {
        for (int i = 0; i < m; i++) {
          const int prev = (i - 1 + m) % m;
          if (!isOuterVertex(mainRing[i]) && isOuterVertex(mainRing[prev])) {
            std::rotate(mainRing.begin(), mainRing.begin() + i, mainRing.end());
            break;
          }
        }
      }
    }
  }

  // 7. Concatenate every closed loop then every open chain into one merged
  // graph, tagging each point's width and inner(0)/outer(1) track.
  zPointArray mergedPositions;
  zIntArray mergedEdges;
  zFloatArray mergedWidths;
  zFloatArray mergedTrackIds;

  auto appendSeq = [&](const std::vector<int> &seq, bool isLoop) {
    const int base = static_cast<int>(mergedPositions.size());
    for (int id : seq) {
      mergedPositions.push_back(rawPositions[id]);
      mergedWidths.push_back(id < int(rawWidths.size()) ? rawWidths[id]
                                                        : 0.04f);
      mergedTrackIds.push_back(isOuterVertex(id) ? 1.0f : 0.0f);
    }
    for (size_t i = 0; i + 1 < seq.size(); i++) {
      mergedEdges.push_back(base + int(i));
      mergedEdges.push_back(base + int(i) + 1);
    }
    if (isLoop && seq.size() >= 3) {
      mergedEdges.push_back(base + int(seq.size()) - 1);
      mergedEdges.push_back(base);
    }
  };

  for (auto &loop : closedLoops)
    appendSeq(loop, true);
  for (auto &chain : openChains)
    appendSeq(chain, false);

  if (mergedPositions.size() < 3 || mergedEdges.size() < 4)
    return;

  PrintPath p;
  zFnGraph fnP(p.graph);
  fnP.create(mergedPositions, mergedEdges);
  p.widths = mergedWidths;
  p.trackIds = mergedTrackIds;
  outPaths.push_back(p);
}

// Fixed-spacing arc-length resample of ioPath.graph, carrying
// widths/heights/normals along. The sampling phase resets at each
// preserved feature/corner vertex (and at the open path's true endpoint)
// so those land exactly on an output point instead of being smoothed
// over; everywhere else, output points are placed at a uniform
// sampleLength apart regardless of the input's original vertex spacing
// (unlike a simple "subdivide long segments" pass, this can also make
// an already-dense input coarser).
bool _resamplePrintPath(PrintPath &ioPath, float sampleLength,
                        float featureAngleTol) {
  zFnGraph fnGraph(ioPath.graph);
  const int N = fnGraph.numVertices();
  if (N < 2 || sampleLength <= 1e-6f)
    return false;

  const bool closed = isClosedGraph(ioPath.graph);

  zPointArray positions;
  fnGraph.getVertexPositions(positions);

  std::vector<int> sequence(N);
  for (int i = 0; i < N; i++)
    sequence[i] = i;

  std::vector<int> featureFlags;
  _detectFeatureAngles(positions, sequence, featureAngleTol, featureFlags);

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

  auto appendPoint = [&](const zPoint &p, float w, float h, const zVector &nrm,
                         float trk) {
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
    return (hasHeights && id < static_cast<int>(heights.size())) ? heights[id]
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
        const float interpTrk =
            (trk0 + t * (trk1 - trk0)) >= 0.5f ? 1.0f : 0.0f;
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
void _filterGlitches(PrintPath &ioPath) {
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

// Orients one layer's ring path in place - CCW if ccw is true, CW if false
// (extraction already produces at most one ring per layer, already
// inner-before-outer ordered).
bool _orientLayer(std::vector<PrintPath> &ringPaths, bool ccw) {
  if (ringPaths.empty())
    return false;

  PrintPath &path = ringPaths[0];
  zFnGraph fnGraph(path.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  if (positions.size() < 3)
    return false;
  const bool isCCW = signedArea2D(positions) >= 0.0;
  if (isCCW != ccw) {
    _reversePrintPath(path);
  }
  return true;
}

} // anonymous namespace

//--------------------------
//---- PUBLIC METHODS
//--------------------------

// Step 1: extracts each layer's isocontour into a single oriented,
// inner-before-outer PrintPath per layer. Parallel across layers
// internally.
ZSPACE_INLINE bool zPrintSyn::extractPrintPaths(
    zObjMeshScalarField &scalarField,
    const std::vector<zScalarArray> &pathFields,
    const std::vector<zScalarArray> &widthFields,
    const std::vector<zScalarArray> &trackFields,
    std::vector<std::vector<PrintPath>> &outPrintPaths2D) {
  const int totalLayers = int(pathFields.size());
  if (totalLayers == 0 || widthFields.size() != size_t(totalLayers) ||
      trackFields.size() != size_t(totalLayers)) {
    m_Status = "failed: mismatched or empty field arrays.";
    printError("[zPrintSyn::extractPrintPaths] " + m_Status);
    return false;
  }

  outPrintPaths2D.assign(totalLayers, std::vector<PrintPath>());

  const auto wallStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel
  {
    zObjMeshScalarField localField = scalarField;

#pragma omp for schedule(dynamic)
    for (int l = 0; l < totalLayers; l++) {
      zPointArray rawPositions;
      zIntArray rawEdges;
      zFloatArray rawWidths;
      zFloatArray rawTrackIds;
      _extractRawIsoContour(localField, pathFields[l], widthFields[l],
                            trackFields[l], rawPositions, rawEdges, rawWidths,
                            rawTrackIds);
      _extractTrackPaths(rawPositions, rawEdges, rawWidths, rawTrackIds, 0.08,
                         1.0, outPrintPaths2D[l]);
    }
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Extracted track paths for " << totalLayers << " layers in "
      << std::fixed << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  printSuccess("[zPrintSyn::extractPrintPaths] " + m_Status);
  return true;
}

// Step 2: runs _resamplePrintPath -> _filterGlitches on every layer's
// paths (in parallel per layer), then forces each layer's ring CCW via
// _orientLayer.
ZSPACE_INLINE bool
zPrintSyn::resamplePrintPath(std::vector<std::vector<PrintPath>> &ioPaths) {
  const int totalLayers = int(ioPaths.size());
  if (totalLayers == 0) {
    m_Status = "failed: No print paths available.";
    printError("[zPrintSyn::resamplePrintPath] " + m_Status);
    return false;
  }

  const float sampleLength =
      (i_Settings.ResampleLength > 1e-4f) ? i_Settings.ResampleLength : 0.005f;

  const auto wallStart = std::chrono::high_resolution_clock::now();
  int totalPaths = 0;

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    for (auto &path : ioPaths[l]) {
      if (_resamplePrintPath(path, sampleLength, i_Settings.FeatureAngleTol)) {
        _filterGlitches(path);
      }
    }
    _orientLayer(ioPaths[l], true);
#pragma omp atomic
    totalPaths += int(ioPaths[l].size());
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Resampled " << totalPaths << " print paths in " << std::fixed
      << std::setprecision(2) << wallTotalMs
      << " ms (sampleLen=" << sampleLength << "mm)";
  m_Status = out.str();
  printSuccess("[zPrintSyn::resamplePrintPath] " + m_Status);
  return true;
}

// Step 3: aligns every layer's ring start point to the same barycentric
// point on the underlying mesh, constrained to the inner track. Seeds from
// the first usable layer by walking i_Settings.SeamOffset from its
// extracted start along the inner segment (clamped so it can't cross into
// outer), parameterizes that point against the seed layer's 2D section
// mesh, then realizes the same parametric point on every other layer's
// mesh and rotates that layer's ring to start at its nearest inner-segment
// point - no further offset there, just the nearest match - tracking the
// same physical spot on the model rather than a fixed 2D region, so it
// stays close layer-to-layer even through taper/twist. Finally cuts every
// ring open at its (now aligned) start via
// _cutGapBeforeStart(i_Settings.SeamSize), so each layer ends as an open
// chain ready for linkPrintPath.
ZSPACE_INLINE bool
zPrintSyn::alignStartPoints(std::vector<std::vector<PrintPath>> &ioPaths,
                            const zObjMeshArray &sectionMeshes2D) {
  const double seamOffset = i_Settings.SeamOffset;
  const double seamSize = i_Settings.SeamSize;
  const int totalLayers = int(ioPaths.size());
  if (totalLayers == 0 || int(sectionMeshes2D.size()) < totalLayers) {
    m_Status = "failed: No print paths or mismatched section meshes.";
    printError("[zPrintSyn::alignStartPoints] " + m_Status);
    return false;
  }

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

  // 1. Seed: walk seamOffset from the current start along the inner
  // segment (index 0 is always inner, as a contiguous prefix, per
  // _extractTrackPaths), clamped so it can't cross into outer.
  PrintPath &seedPath = ioPaths[seedLayer][0];
  {
    zFnGraph fnGraph(seedPath.graph);
    zPointArray positions;
    fnGraph.getVertexPositions(positions);
    const std::vector<double> cumLen = _cumulativeArcLength(positions);
    const double innerLen =
        _innerSegmentLength(positions, seedPath.trackIds, cumLen);
    const double clampedWalk =
        (std::max)(0.0,
                   (std::min)(seamOffset, (std::max)(0.0, innerLen - 1e-4)));
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
    seedPts.push_back(
        positions[0]); // duplicate: parameterizeGraph needs an edge
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
    if (!realizeParameterGraph(
            seedParams, const_cast<zObjMesh &>(sectionMeshes2D[l]), realized))
      continue;
    zFnGraph fnRealized(realized);
    zPointArray realizedPts;
    fnRealized.getVertexPositions(realizedPts);
    if (realizedPts.empty())
      continue;

    double arcPos = 0.0;
    if (!_nearestInnerArcPos(path, realizedPts[0], arcPos))
      continue;

    if (_insertAndRotateToStart(path, arcPos)) {
#pragma omp atomic
      alignedCount++;
    }
  }

  // 3. Cut every layer's ring open right before its (now aligned) start
  // point, leaving a seamSize gap between the new end point and the
  // start - ready to be chained layer-to-layer. Independent per layer
  // (each only touches its own ioPaths[l][0]), so this runs in parallel
  // too, same as step 2.
#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    if (ioPaths[l].empty())
      continue;
    _cutGapBeforeStart(ioPaths[l][0], seamSize);
  }

  if (alignedCount == 0) {
    m_Status = "failed: No layers could be aligned.";
    printError("[zPrintSyn::alignStartPoints] " + m_Status);
    return false;
  }

  std::ostringstream out;
  out << "Aligned start points for " << alignedCount << " of " << totalLayers
      << " layers";
  m_Status = out.str();
  printSuccess("[zPrintSyn::alignStartPoints] " + m_Status);
  return true;
}

// Step 5: concatenates every layer's single (post-alignStartPoints, open)
// PrintPath end-to-start into one continuous bottom-to-top PrintPath.
// Simple sequential concatenation already connects each layer's last
// point to the next layer's first point via the edge between them.
// Sequential (order-dependent), no internal parallelism.
ZSPACE_INLINE bool zPrintSyn::linkPrintPath(
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
                                               // alignStartPoints

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
        combinedHeights.insert(combinedHeights.end(), layerPath.heights.begin(),
                               layerPath.heights.end());
      else
        combinedHeights.insert(combinedHeights.end(), positions.size(), 0.0f);
    }
    if (hasNormals) {
      if (!layerPath.normals.empty())
        combinedNormals.insert(combinedNormals.end(), layerPath.normals.begin(),
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
        combinedTrackIds.insert(combinedTrackIds.end(), positions.size(), 0.0f);
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

  std::ostringstream out;
  out << "Linked " << layeredPaths.size() << " layers into one print path ("
      << total << " points)";
  m_Status = out.str();
  printSuccess("[zPrintSyn::linkPrintPath] " + m_Status);
  return true;
}

// Step 4: projects 2D PrintPaths back onto the 3D section meshes.
// Parallel across layers internally.
ZSPACE_INLINE bool zPrintSyn::projectPrintPath(
    const std::vector<std::vector<PrintPath>> &in2DPaths,
    const zObjMeshArray &sectionMeshes2D, const zObjMeshArray &sectionMeshes,
    const std::vector<zScalarArray> &vScalar, float layerSpacing,
    std::vector<std::vector<PrintPath>> &out3DPaths) {
  const int totalLayers = int(in2DPaths.size());
  if (totalLayers == 0 || sectionMeshes2D.empty() || sectionMeshes.empty()) {
    m_Status = "failed: No 2D print paths or section meshes available for "
              "projection.";
    printError("[zPrintSyn::projectPrintPath] " + m_Status);
    return false;
  }

  out3DPaths.assign(totalLayers, std::vector<PrintPath>());
  for (int l = 0; l < totalLayers; l++) {
    out3DPaths[l].assign(in2DPaths[l].size(), PrintPath());
  }

  zFloatArray meshVertexHeights;
  if (!vScalar.empty() && totalLayers > 1) {
    meshVertexHeights.resize(vScalar.size());
    const float denom = float(totalLayers - 1);
    for (size_t k = 0; k < vScalar.size(); k++) {
      meshVertexHeights[k] = (!vScalar[k].empty())
                                 ? float(vScalar[k].back() / denom)
                                 : layerSpacing;
    }
  }

  int totalPaths = 0;
  const auto wallStart = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic)
  for (int l = 0; l < totalLayers; l++) {
    for (size_t r = 0; r < in2DPaths[l].size(); r++) {
      out3DPaths[l][r].graph = in2DPaths[l][r].graph;
      out3DPaths[l][r].widths = in2DPaths[l][r].widths;
      out3DPaths[l][r].trackIds = in2DPaths[l][r].trackIds;

      zFnGraph fnContour(out3DPaths[l][r].graph);
      if (fnContour.numVertices() > 0 && l < int(sectionMeshes2D.size()) &&
          l < int(sectionMeshes.size())) {
        remapGraph(
            out3DPaths[l][r].graph, const_cast<zObjMesh &>(sectionMeshes2D[l]),
            const_cast<zObjMesh &>(sectionMeshes[l]), &out3DPaths[l][r].normals,
            &meshVertexHeights, &out3DPaths[l][r].heights);
      }
    }
#pragma omp atomic
    totalPaths += int(in2DPaths[l].size());
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Projected " << totalPaths << " 3D PrintPaths in " << std::fixed
      << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  printSuccess("[zPrintSyn::projectPrintPath] " + m_Status);
  return true;
}

} // namespace zSpace
