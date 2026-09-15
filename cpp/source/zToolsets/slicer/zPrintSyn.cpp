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

// Orients one layer's 1-2 ring paths into a single closed loop.
ZSPACE_INLINE bool
zPrintSyn::orientAndWeldLayer(std::vector<PrintPath> &ringPaths) {
  if (ringPaths.empty())
    return false;

  // Single ring (single-perimeter boundary): just force it CCW.
  if (ringPaths.size() == 1) {
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

  // Two rings: locate inner/outer explicitly (vector order isn't
  // guaranteed to be [inner, outer]).
  int innerIdx = -1, outerIdx = -1;
  for (size_t i = 0; i < ringPaths.size(); i++) {
    if (ringPaths[i].isOuterTrack)
      outerIdx = int(i);
    else
      innerIdx = int(i);
  }
  if (innerIdx < 0 || outerIdx < 0)
    return false; // not a clean inner+outer pair - leave as-is

  PrintPath &inner = ringPaths[innerIdx];
  PrintPath &outer = ringPaths[outerIdx];

  zFnGraph fnInner(inner.graph);
  zPointArray innerPos;
  fnInner.getVertexPositions(innerPos);
  if (innerPos.size() < 2)
    return false;

  // Force inner counter-clockwise.
  if (_signedArea2D(innerPos) < 0.0) {
    _reversePrintPath(inner);
    fnInner.getVertexPositions(innerPos);
  }

  zFnGraph fnOuter(outer.graph);
  zPointArray outerPos;
  fnOuter.getVertexPositions(outerPos);
  if (outerPos.size() < 2)
    return false;

  // Orient outer so it starts at the end where inner left off (the two
  // chains were already made to touch there during extraction).
  const zPoint &innerEnd = innerPos.back();
  const double distToFront = vecDistSq(innerEnd, outerPos.front());
  const double distToBack = vecDistSq(innerEnd, outerPos.back());
  if (distToBack < distToFront) {
    _reversePrintPath(outer);
    fnOuter.getVertexPositions(outerPos);
  }

  // Weld: inner's points first, then outer's, closed back to inner's start.
  const bool hasHeights = !inner.heights.empty() || !outer.heights.empty();
  const bool hasNormals = !inner.normals.empty() || !outer.normals.empty();

  zPointArray combinedPositions;
  zFloatArray combinedWidths;
  zFloatArray combinedHeights;
  zVectorArray combinedNormals;

  combinedPositions.insert(combinedPositions.end(), innerPos.begin(),
                           innerPos.end());
  combinedPositions.insert(combinedPositions.end(), outerPos.begin(),
                           outerPos.end());
  combinedWidths.insert(combinedWidths.end(), inner.widths.begin(),
                        inner.widths.end());
  combinedWidths.insert(combinedWidths.end(), outer.widths.begin(),
                        outer.widths.end());
  if (hasHeights) {
    combinedHeights.insert(combinedHeights.end(), inner.heights.begin(),
                           inner.heights.end());
    combinedHeights.insert(combinedHeights.end(), outer.heights.begin(),
                           outer.heights.end());
  }
  if (hasNormals) {
    combinedNormals.insert(combinedNormals.end(), inner.normals.begin(),
                           inner.normals.end());
    combinedNormals.insert(combinedNormals.end(), outer.normals.begin(),
                           outer.normals.end());
  }

  const int total = static_cast<int>(combinedPositions.size());
  if (total < 3)
    return false;

  zIntArray combinedEdges;
  combinedEdges.reserve(total * 2);
  for (int i = 0; i + 1 < total; i++) {
    combinedEdges.push_back(i);
    combinedEdges.push_back(i + 1);
  }
  combinedEdges.push_back(total - 1);
  combinedEdges.push_back(0);

  PrintPath merged;
  zFnGraph fnMerged(merged.graph);
  fnMerged.create(combinedPositions, combinedEdges);
  merged.widths = combinedWidths;
  if (hasHeights)
    merged.heights = combinedHeights;
  if (hasNormals)
    merged.normals = combinedNormals;
  merged.isOuterTrack = false; // no longer a single track once merged

  ringPaths.clear();
  ringPaths.push_back(merged);
  return true;
}

// Cuts a gap of width `width` out of ioPath's closed graph, centered on the
// point of ioPath nearest to nearPoint.
ZSPACE_INLINE bool zPrintSyn::breakPrintPath(PrintPath &ioPath,
                                             const zPoint &nearPoint,
                                             float width) {
  if (width <= 1e-6f)
    return false;
  if (!isClosedGraph(ioPath.graph))
    return false; // only meaningful for a closed loop

  zFnGraph fnGraph(ioPath.graph);
  zPointArray positions;
  fnGraph.getVertexPositions(positions);
  const int numV = static_cast<int>(positions.size());
  if (numV < 3)
    return false;

  // 1. Cumulative arc length around the loop (segment i: point i -> point
  // (i+1)%numV, including the closing segment).
  std::vector<double> cumLen(numV + 1, 0.0);
  for (int i = 0; i < numV; i++) {
    cumLen[i + 1] =
        cumLen[i] + vecSub(positions[(i + 1) % numV], positions[i]).length();
  }
  const double L = cumLen[numV];
  if (L <= 1e-9 || double(width) >= L)
    return false;

  // 2. Nearest point on the loop to nearPoint -> its global arc-length
  // position.
  double bestDistSq = (std::numeric_limits<double>::max)();
  double nearestArcPos = 0.0;
  for (int i = 0; i < numV; i++) {
    const zPoint &p0 = positions[i];
    const zPoint &p1 = positions[(i + 1) % numV];
    zVector seg = vecSub(p1, p0);
    const double segLenSq = vecDot(seg, seg);
    double t = 0.0;
    if (segLenSq > 1e-12) {
      t = vecDot(vecSub(nearPoint, p0), seg) / segLenSq;
      t = (std::max)(0.0, (std::min)(1.0, t));
    }
    zPoint proj(p0.x + t * seg.x, p0.y + t * seg.y, p0.z + t * seg.z);
    const double distSq = vecDistSq(nearPoint, proj);
    if (distSq < bestDistSq) {
      bestDistSq = distSq;
      nearestArcPos = cumLen[i] + t * (cumLen[i + 1] - cumLen[i]);
    }
  }

  // 3. Forward/backward cut points, 0.5*width along the loop from the
  // nearest point (wraps around since the loop is closed).
  const double half = 0.5 * double(width);
  double aPos = std::fmod(nearestArcPos + half, L);
  if (aPos < 0.0)
    aPos += L;
  double bPos = std::fmod(nearestArcPos - half + L, L);
  if (bPos < 0.0)
    bPos += L;

  const bool hasHeights = !ioPath.heights.empty();
  const bool hasNormals = !ioPath.normals.empty();

  // Interpolates position/width/height/normal at global arc-length s.
  auto sampleAt = [&](double s, zPoint &outP, float &outW, float &outH,
                      zVector &outN) {
    int seg = 0;
    for (; seg < numV; seg++) {
      if (s <= cumLen[seg + 1] || seg == numV - 1)
        break;
    }
    const int i0 = seg, i1 = (seg + 1) % numV;
    const double segStart = cumLen[seg], segEnd = cumLen[seg + 1];
    const double segLen = segEnd - segStart;
    const double t = (segLen > 1e-9)
                         ? (std::max)(0.0, (std::min)(1.0,
                                                      (s - segStart) / segLen))
                         : 0.0;
    const zPoint &p0 = positions[i0];
    const zPoint &p1 = positions[i1];
    outP = zPoint(p0.x + t * (p1.x - p0.x), p0.y + t * (p1.y - p0.y),
                 p0.z + t * (p1.z - p0.z));
    const float w0 =
        (i0 < int(ioPath.widths.size())) ? ioPath.widths[i0] : 0.04f;
    const float w1 =
        (i1 < int(ioPath.widths.size())) ? ioPath.widths[i1] : 0.04f;
    outW = static_cast<float>(w0 + t * (w1 - w0));
    if (hasHeights) {
      const float h0 =
          (i0 < int(ioPath.heights.size())) ? ioPath.heights[i0] : 0.0f;
      const float h1 =
          (i1 < int(ioPath.heights.size())) ? ioPath.heights[i1] : 0.0f;
      outH = static_cast<float>(h0 + t * (h1 - h0));
    } else {
      outH = 0.0f;
    }
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
      outN = n;
    } else {
      outN = zVector(0, 0, 1);
    }
  };

  zPointArray newPositions;
  zFloatArray newWidths;
  zFloatArray newHeights;
  zVectorArray newNormals;

  zPoint pA;
  float wA, hA;
  zVector nA;
  sampleAt(aPos, pA, wA, hA, nA);
  newPositions.push_back(pA);
  newWidths.push_back(wA);
  if (hasHeights)
    newHeights.push_back(hA);
  if (hasNormals)
    newNormals.push_back(nA);

  // Keep every original vertex strictly between A and B going the long way
  // around (i.e. NOT through the gap that contains the original nearest
  // point).
  const double keepSpan = std::fmod(bPos - aPos + L, L);
  for (int i = 0; i < numV; i++) {
    const double relPos = std::fmod(cumLen[i] - aPos + L, L);
    if (relPos > 1e-9 && relPos < keepSpan - 1e-9) {
      newPositions.push_back(positions[i]);
      newWidths.push_back((i < int(ioPath.widths.size())) ? ioPath.widths[i]
                                                           : 0.04f);
      if (hasHeights)
        newHeights.push_back(
            (i < int(ioPath.heights.size())) ? ioPath.heights[i] : 0.0f);
      if (hasNormals)
        newNormals.push_back((i < int(ioPath.normals.size()))
                                 ? ioPath.normals[i]
                                 : zVector(0, 0, 1));
    }
  }

  zPoint pB;
  float wB, hB;
  zVector nB;
  sampleAt(bPos, pB, wB, hB, nB);
  newPositions.push_back(pB);
  newWidths.push_back(wB);
  if (hasHeights)
    newHeights.push_back(hB);
  if (hasNormals)
    newNormals.push_back(nB);

  if (newPositions.size() < 2)
    return false;

  const int newN = static_cast<int>(newPositions.size());
  zIntArray newEdges;
  newEdges.reserve((newN - 1) * 2);
  for (int i = 0; i + 1 < newN; i++) {
    newEdges.push_back(i);
    newEdges.push_back(i + 1);
  }

  fnGraph.clear();
  fnGraph.create(newPositions, newEdges);
  ioPath.widths = newWidths;
  if (hasHeights)
    ioPath.heights = newHeights;
  if (hasNormals)
    ioPath.normals = newNormals;
  return true;
}

// Concatenates every layer's single PrintPath end-to-start into one
// continuous bottom-to-top PrintPath.
ZSPACE_INLINE bool zPrintSyn::chainLayersIntoOnePath(
    const std::vector<std::vector<PrintPath>> &layeredPaths,
    PrintPath &outCombined) {
  outCombined = PrintPath();

  bool hasHeights = false, hasNormals = false;
  for (const auto &ringPaths : layeredPaths) {
    if (!ringPaths.empty()) {
      if (!ringPaths[0].heights.empty())
        hasHeights = true;
      if (!ringPaths[0].normals.empty())
        hasNormals = true;
    }
  }

  zPointArray combinedPositions;
  zFloatArray combinedWidths;
  zFloatArray combinedHeights;
  zVectorArray combinedNormals;

  for (const auto &ringPaths : layeredPaths) {
    if (ringPaths.empty())
      continue;
    const PrintPath &layerPath = ringPaths[0]; // one path per layer, post
                                               // orientAndWeldLayer

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
  return true;
}

} // namespace zSpace
