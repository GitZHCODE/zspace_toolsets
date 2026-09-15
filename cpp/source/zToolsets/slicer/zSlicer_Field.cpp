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

#include <headers/zToolsets/slicer/zSlicer.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
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

// Builds a perpendicular trim slot segment at ratio t along a boundary/rib
// graph.
bool _createPerpendicularSlot(const zObjGraph &sourceGraph, zObjGraph &outGraph,
                              zScalar t, zScalar trimLength) {
  zFnGraph fnSource(const_cast<zObjGraph &>(sourceGraph));
  zPointArray positions;
  zIntArray edgeConnects;
  fnSource.getVertexPositions(positions);
  fnSource.getEdgeData(edgeConnects);

  zFnGraph(outGraph).clear();
  if (positions.size() < 2 || edgeConnects.size() < 2)
    return false;

  const int numEdges = static_cast<int>(edgeConnects.size()) / 2;
  const zScalar ratio = (std::max)(zScalar(0.0), (std::min)(zScalar(1.0), t));

  // 1. Construct lightweight adjacency list
  std::vector<std::vector<std::pair<int, int>>> adj(positions.size());
  for (int e = 0; e < numEdges; e++) {
    int u = edgeConnects[e * 2], v = edgeConnects[e * 2 + 1];
    if (u >= 0 && v >= 0 && u < (int)positions.size() &&
        v < (int)positions.size() &&
        vecDistSq(positions[u], positions[v]) > 1e-12) {
      adj[u].push_back({v, e});
      adj[v].push_back({u, e});
    }
  }

  std::vector<bool> visitedEdge(numEdges, false),
      visitedVert(positions.size(), false);
  zPointArray trimPositions;
  zIntArray trimEdges;

  // 2. Traverse connected components and generate perpendicular slot at ratio
  for (size_t seed = 0; seed < positions.size(); seed++) {
    if (visitedVert[seed] || adj[seed].empty())
      continue;

    // Collect connected component vertices
    std::vector<int> comp, q = {int(seed)};
    visitedVert[seed] = true;
    while (!q.empty()) {
      int curr = q.back();
      q.pop_back();
      comp.push_back(curr);
      for (const auto &nb : adj[curr]) {
        if (!visitedVert[nb.first]) {
          visitedVert[nb.first] = true;
          q.push_back(nb.first);
        }
      }
    }

    // Find endpoint (vertex of degree 1; if closed loop, start at first vertex)
    int start = comp[0];
    for (int v : comp) {
      if (adj[v].size() == 1) {
        start = v;
        break;
      }
    }

    // Walk polyline chain and build ordered positions and cumulative arc
    // lengths
    std::vector<zPoint> chainPos = {positions[start]};
    std::vector<zScalar> chainScalars = {0.0};
    double totalLen = 0.0;
    int cur = start, prev = -1;

    while (true) {
      int next = -1, edgeId = -1;
      for (const auto &nb : adj[cur]) {
        if (nb.first != prev && !visitedEdge[nb.second]) {
          next = nb.first;
          edgeId = nb.second;
          break;
        }
      }
      if (next == -1)
        break;
      visitedEdge[edgeId] = true;
      totalLen += std::sqrt(vecDistSq(positions[cur], positions[next]));
      chainPos.push_back(positions[next]);
      chainScalars.push_back(static_cast<zScalar>(totalLen));
      prev = cur;
      cur = next;
    }

    if (totalLen <= 1e-6 || chainPos.size() < 2)
      continue;

    // Locate midpoint and perpendicular direction using shared samplePolyline
    zPoint mid;
    zVector dir;
    if (samplePolyline(chainPos.data(), chainScalars.data(), chainPos.size(),
                       ratio, mid, &dir)) {
      zVector perp(-dir.y, dir.x, 0.0);
      int baseId = static_cast<int>(trimPositions.size());
      trimPositions.push_back(mid + perp * trimLength);
      trimPositions.push_back(mid - perp * trimLength);
      trimEdges.push_back(baseId);
      trimEdges.push_back(baseId + 1);
    }
  }

  if (trimPositions.empty())
    return false;
  zFnGraph(outGraph).create(trimPositions, trimEdges);
  return true;
}

struct FastSegment {
  float x0, y0;
  float dx, dy;
  float lenSqInv;
  float minX, maxX;
  float minY, maxY;
};

// Converts a graph's edges into flat FastSegment records for fast distance
// queries.
inline void computeFastSegments(const zObjGraph &graph,
                                std::vector<FastSegment> &outSegments) {
  outSegments.clear();
  zFnGraph fn(const_cast<zObjGraph &>(graph));
  zPointArray positions;
  zIntArray edges;
  fn.getVertexPositions(positions);
  fn.getEdgeData(edges);
  if (positions.empty() || edges.empty())
    return;

  outSegments.reserve(edges.size() / 2);
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    int v0 = edges[e];
    int v1 = edges[e + 1];
    if (v0 < 0 || v0 >= (int)positions.size() || v1 < 0 ||
        v1 >= (int)positions.size())
      continue;

    float x0 = static_cast<float>(positions[v0].x);
    float y0 = static_cast<float>(positions[v0].y);
    float x1 = static_cast<float>(positions[v1].x);
    float y1 = static_cast<float>(positions[v1].y);

    float dx = x1 - x0;
    float dy = y1 - y0;
    float lenSq = dx * dx + dy * dy;
    if (lenSq < 1e-12f)
      continue;

    FastSegment seg;
    seg.x0 = x0;
    seg.y0 = y0;
    seg.dx = dx;
    seg.dy = dy;
    seg.lenSqInv = 1.0f / lenSq;
    seg.minX = (std::min)(x0, x1);
    seg.maxX = (std::max)(x0, x1);
    seg.minY = (std::min)(y0, y1);
    seg.maxY = (std::max)(y0, y1);
    outSegments.push_back(seg);
  }
}

// Computes per-point unsigned distance to the nearest FastSegment, minus
// offset.
inline void computeFastDistanceField(const zPoint *points, int numPoints,
                                     const std::vector<FastSegment> &segments,
                                     zScalarArray &outDistances,
                                     float offset = 0.0f) {
  outDistances.resize(numPoints);
  if (segments.empty()) {
    std::fill(outDistances.begin(), outDistances.end(),
              (std::numeric_limits<zScalar>::max)());
    return;
  }

  const int numSegs = static_cast<int>(segments.size());

  for (int i = 0; i < numPoints; i++) {
    const float px = static_cast<float>(points[i].x);
    const float py = static_cast<float>(points[i].y);
    float minD2 = (std::numeric_limits<float>::max)();

    for (int k = 0; k < numSegs; k++) {
      const FastSegment &seg = segments[k];

      float dboxX = (std::max)(0.0f, (std::max)(seg.minX - px, px - seg.maxX));
      float dboxY = (std::max)(0.0f, (std::max)(seg.minY - py, py - seg.maxY));
      if (dboxX * dboxX + dboxY * dboxY >= minD2)
        continue;

      float t =
          ((px - seg.x0) * seg.dx + (py - seg.y0) * seg.dy) * seg.lenSqInv;
      t = (std::max)(0.0f, (std::min)(1.0f, t));
      float projX = seg.x0 + t * seg.dx;
      float projY = seg.y0 + t * seg.dy;
      float dist2 = (px - projX) * (px - projX) + (py - projY) * (py - projY);
      if (dist2 < minD2)
        minD2 = dist2;
    }

    outDistances[i] = static_cast<zScalar>(std::sqrt(minD2) - offset);
  }
}

// Computes per-point inside/outside sign (-1/1) via winding number against
// polygonGraph.
inline void computeFastPolygon(const zPoint *points, int numPoints,
                               const zObjGraph &polygonGraph,
                               zScalarArray &outPolyField) {
  outPolyField.resize(numPoints);
  if (numPoints == 0)
    return;

  zFnGraph fn(const_cast<zObjGraph &>(polygonGraph));
  zPointArray positions;
  zIntArray edges;
  fn.getVertexPositions(positions);
  fn.getEdgeData(edges);

  if (positions.empty() || edges.empty()) {
    std::fill(outPolyField.begin(), outPolyField.end(), 1.0);
    return;
  }

  struct PolyEdge {
    float x0, y0;
    float x1, y1;
    float minY, maxY;
  };

  std::vector<PolyEdge> polyEdges;
  polyEdges.reserve(edges.size() / 2);

  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    int v0 = edges[e];
    int v1 = edges[e + 1];
    if (v0 < 0 || v0 >= (int)positions.size() || v1 < 0 ||
        v1 >= (int)positions.size())
      continue;

    PolyEdge pe;
    pe.x0 = static_cast<float>(positions[v0].x);
    pe.y0 = static_cast<float>(positions[v0].y);
    pe.x1 = static_cast<float>(positions[v1].x);
    pe.y1 = static_cast<float>(positions[v1].y);
    pe.minY = (std::min)(pe.y0, pe.y1);
    pe.maxY = (std::max)(pe.y0, pe.y1);
    polyEdges.push_back(pe);
  }

  const int numEdges = static_cast<int>(polyEdges.size());

  for (int i = 0; i < numPoints; i++) {
    const float px = static_cast<float>(points[i].x);
    const float py = static_cast<float>(points[i].y);

    int windingNumber = 0;

    for (int e = 0; e < numEdges; e++) {
      const PolyEdge &edge = polyEdges[e];

      if (py < edge.minY || py > edge.maxY)
        continue;

      if (edge.y0 <= py) {
        if (edge.y1 > py) {
          float cross = (edge.x1 - edge.x0) * (py - edge.y0) -
                        (px - edge.x0) * (edge.y1 - edge.y0);
          if (cross > 0.0f)
            windingNumber++;
        }
      } else {
        if (edge.y1 <= py) {
          float cross = (edge.x1 - edge.x0) * (py - edge.y0) -
                        (px - edge.x0) * (edge.y1 - edge.y0);
          if (cross < 0.0f)
            windingNumber--;
        }
      }
    }

    outPolyField[i] = (windingNumber != 0) ? -1.0 : 1.0;
  }
}

// Groups a boundary/bracing graph with its print width, offset and trim ratio.
struct SegmentGraph {
  zObjGraph graph;
  int type = 0;         // 0: Outer Boundary, 1: Internal Bracing
  zScalar width = 0.0;  // Print bead width
  zScalar offset = 0.0; // Distance offset
  zScalar trim = 0.0;   // Trim slot ratio

  SegmentGraph() = default;
  SegmentGraph(const zObjGraph &inGraph, int inType, zScalar inWidth,
               zScalar inOffset, zScalar inTrim = 0.0)
      : graph(inGraph), type(inType), width(inWidth), offset(inOffset),
        trim(inTrim) {}
};

// Computes per-field-vertex distance to each segment graph.
void _computeSegmentField(zObjMeshScalarField &field,
                          const std::vector<SegmentGraph> &segmentGraphs,
                          std::vector<zScalarArray> &outSegDistances) {
  zFnMesh fnMesh(field);
  const zPoint *points = fnMesh.getRawVertexPositions();
  const int n = fnMesh.numVertices();
  const int numSegs = static_cast<int>(segmentGraphs.size());
  outSegDistances.assign(
      numSegs, zScalarArray(n, (std::numeric_limits<zScalar>::max)()));
  if (!points || n <= 0)
    return;

  std::vector<FastSegment> fastEdges;
  for (int k = 0; k < numSegs; k++) {
    computeFastSegments(segmentGraphs[k].graph, fastEdges);
    computeFastDistanceField(points, n, fastEdges, outSegDistances[k]);
  }
}

// Assigns each field vertex the print width of its dominating (nearest)
// segment.
void _computeWidthField(zObjMeshScalarField &field,
                        const std::vector<zScalarArray> &segDistances,
                        const std::vector<SegmentGraph> &segmentGraphs,
                        bool doublelayerboundary, zScalar overlapWidth,
                        zScalarArray *outWidthField = nullptr,
                        zScalarArray *outIdField = nullptr) {
  zFnMeshScalarField fnField(field);
  const int n = fnField.numFieldValues();
  zScalarArray widthValues(n, 0.0);
  zScalarArray idValues(n, 0.0);

  if (segDistances.empty() || segmentGraphs.empty() ||
      segDistances.size() != segmentGraphs.size()) {
    return;
  }

  const int numSegs = static_cast<int>(segmentGraphs.size());
  std::vector<zScalar> actualPrintWidths(numSegs);
  for (int k = 0; k < numSegs; k++) {
    const bool isOuter = (segmentGraphs[k].type == 0);
    if (isOuter && !doublelayerboundary) {
      actualPrintWidths[k] = segmentGraphs[k].width;
    } else {
      actualPrintWidths[k] =
          (std::max)(zScalar(1e-4), zScalar(0.5 * segmentGraphs[k].width +
                                            0.5 * overlapWidth));
    }
  }

  std::vector<zScalar> minMetrics(n, (std::numeric_limits<zScalar>::max)());

  for (int k = 0; k < numSegs; k++) {
    if (static_cast<int>(segDistances[k].size()) < n)
      continue;

    const zScalar w =
        (segmentGraphs[k].width > 1e-4) ? segmentGraphs[k].width : 1e-4;
    const zScalar off = segmentGraphs[k].offset;
    const zScalarArray &tempDist = segDistances[k];

    for (int i = 0; i < n; i++) {
      // 1. Normalized distance metric based on width and offset
      const zScalar metric = (tempDist[i] - off) / w;
      // 2. Identify dominating segment ID and assign corresponding print width
      if (metric < minMetrics[i]) {
        minMetrics[i] = metric;
        idValues[i] = static_cast<zScalar>(k);
        widthValues[i] = actualPrintWidths[k];
      }
    }
  }

  fnField.setFieldValues(widthValues);
  if (outWidthField)
    *outWidthField = widthValues;
  if (outIdField)
    *outIdField = idValues;
}

// Convenience wrapper: builds FastSegments for graph then computes offset
// distance field.
inline void _evaluateGraphDistance(const zPoint *points, int n,
                                   const zObjGraph &graph, float halfWidth,
                                   std::vector<FastSegment> &fastEdgesBuffer,
                                   zScalarArray &outDist) {
  if (!points || n <= 0)
    return;
  computeFastSegments(graph, fastEdgesBuffer);
  computeFastDistanceField(points, n, fastEdgesBuffer, outDist, halfWidth);
}

// Combines outer wall and internal rib fields (with trim slots) into the final
// SDF path field.
void _computePathField(zObjMeshScalarField &field, const zObjGraph &flatGraph,
                       const std::vector<SegmentGraph> &segmentGraphs,
                       const std::vector<zScalarArray> &segDistances,
                       bool doublelayerboundary, zScalar overlapWidth,
                       zObjGraph *outTrimSlotsGraph = nullptr,
                       zScalarArray *outPathField = nullptr,
                       zScalarArray *outTrackIdField = nullptr) {
  zFnMeshScalarField fnField(field);
  zFnMesh fnMesh(field);
  const zPoint *points = fnMesh.getRawVertexPositions();
  const int n = fnField.numFieldValues();
  zScalarArray pathValues(n, 0.0);
  // Track tag: 1 = outer, 0 = inner. Single-perimeter paths have no
  // inner/outer distinction, so they default to "inner" everywhere.
  zScalarArray trackIdValues(n, 0.0);

  if (!points || n <= 0 || segmentGraphs.empty() || segDistances.empty() ||
      segmentGraphs.size() != segDistances.size()) {
    return;
  }

  const int nS = static_cast<int>(segmentGraphs.size());
  zScalarArray polyField;
  computeFastPolygon(points, n, flatGraph, polyField);

  std::vector<int> outerIndices;
  std::vector<int> ribIndices;
  std::vector<zObjGraph> allSlotGraphs;

  for (int k = 0; k < nS; k++) {
    if (segmentGraphs[k].type == 0)
      outerIndices.push_back(k);
    else
      ribIndices.push_back(k);
  }

  std::vector<FastSegment> fastSlotEdges;

  if (doublelayerboundary) {
    // 1-2. Outer ring field and inner ring field: the same signed
    // polygon-distance quantity per outer segment, offset by that
    // segment's own overlap-adjusted track width, so track 1 (outer)
    // always starts at offset 0 and track 2 (inner)'s far edge always
    // lands at offset w.
    zScalarArray outerRingField(n, 0.0);
    zScalarArray innerRingField(n, 0.0);
    for (int i = 0; i < n; i++) {
      const bool isInside = (polyField[i] <= 0.0);
      double maxOuter = -(std::numeric_limits<double>::max)();
      double maxInner = -(std::numeric_limits<double>::max)();
      for (int k : outerIndices) {
        const zScalar w = segmentGraphs[k].width;
        const zScalar beadW =
            (std::max)(zScalar(1e-4), zScalar(0.5 * w + 0.5 * overlapWidth));
        const zScalar d = segDistances[k][i];
        const zScalar phi1 = (isInside ? -d : d) + 0.5 * beadW;
        const zScalar phi2 = (isInside ? -d : d) + (w - 0.5 * beadW);
        if (phi1 > maxOuter)
          maxOuter = phi1;
        if (phi2 > maxInner)
          maxInner = phi2;
      }
      outerRingField[i] = maxOuter;
      innerRingField[i] = maxInner;
    }

    // 3. Carve internal bracing out of the inner ring, using the exact
    // same combination method the single-perimeter path below uses to
    // combine ribs with a solid field (the inner ring is itself a solid
    // region, same as the single-perimeter outer wall field).
    zScalarArray innerRibsField(n, (std::numeric_limits<zScalar>::max)());
    bool hasValidRibs = false;

    for (int k : ribIndices) {
      const zScalar w = segmentGraphs[k].width;
      const zScalar ribBeadW =
          (std::max)(zScalar(1e-4), zScalar(0.5 * w + 0.5 * overlapWidth));
      const zScalar ribOffset = 0.5 * w - 0.5 * ribBeadW;
      zScalarArray singleRibField(n, 0.0);
      for (int i = 0; i < n; i++) {
        singleRibField[i] = segDistances[k][i] - ribOffset;
      }

      const zScalar trimRatio = segmentGraphs[k].trim;
      zScalarArray singleRibTrimmed = singleRibField;
      if (trimRatio > 1e-4 && trimRatio < 1.0 - 1e-4) {
        zObjGraph slotGraph;
        if (_createPerpendicularSlot(segmentGraphs[k].graph, slotGraph,
                                     trimRatio, w * 1.5)) {
          allSlotGraphs.push_back(slotGraph);
          zScalarArray tempSlotDist;
          const float slotHalfWidth = static_cast<float>(0.5 * ribBeadW);
          _evaluateGraphDistance(points, n, slotGraph, slotHalfWidth,
                                 fastSlotEdges, tempSlotDist);
          fnField.boolean_subtract(singleRibField, tempSlotDist,
                                   singleRibTrimmed, false);
        }
      }

      if (!hasValidRibs) {
        innerRibsField = singleRibTrimmed;
        hasValidRibs = true;
      } else {
        for (int i = 0; i < n; i++) {
          if (singleRibTrimmed[i] < innerRibsField[i]) {
            innerRibsField[i] = singleRibTrimmed[i];
          }
        }
      }
    }

    zScalarArray bracedInnerRing = innerRingField;
    if (hasValidRibs) {
      fnField.boolean_subtract(innerRingField, innerRibsField, bracedInnerRing,
                               false);
    }

    // 4. Outer ring minus (braced) inner ring -> the double-wall band.
    // Tag each grid point by which ring actually dominates there, so the
    // outer/inner distinction survives even though both rings now live on
    // the same combined field.
    zScalarArray ringBand;
    fnField.boolean_subtract(outerRingField, bracedInnerRing, ringBand, false);
    for (int i = 0; i < n; i++) {
      trackIdValues[i] =
          (outerRingField[i] >= -bracedInnerRing[i]) ? 1.0 : 0.0;
    }

    // 5. Cut the boundary trim slot through the combined band, fusing the
    // outer and inner loops into one connected chain. The tag is kept
    // everywhere else; the newly-cut opening itself reads as outer.
    zScalarArray outerSlotDist(n, (std::numeric_limits<zScalar>::max)());
    bool hasOuterSlot = false;

    for (int k : outerIndices) {
      const zScalar trimRatio = segmentGraphs[k].trim;
      if (trimRatio > 1e-4 && trimRatio < 1.0 - 1e-4) {
        zObjGraph slotGraph;
        const zScalar w = segmentGraphs[k].width;
        const zScalar beadW =
            (std::max)(zScalar(1e-4), zScalar(0.5 * w + 0.5 * overlapWidth));
        if (_createPerpendicularSlot(segmentGraphs[k].graph, slotGraph,
                                     trimRatio, w * 1.5)) {
          allSlotGraphs.push_back(slotGraph);
          zScalarArray tempSlotDist;
          const float slotHalfWidth = static_cast<float>(0.5 * beadW);
          _evaluateGraphDistance(points, n, slotGraph, slotHalfWidth,
                                 fastSlotEdges, tempSlotDist);

          if (!hasOuterSlot) {
            outerSlotDist = tempSlotDist;
            hasOuterSlot = true;
          } else {
            for (int i = 0; i < n; i++) {
              if (tempSlotDist[i] < outerSlotDist[i])
                outerSlotDist[i] = tempSlotDist[i];
            }
          }
        }
      }
    }

    if (hasOuterSlot) {
      for (int i = 0; i < n; i++) {
        if (ringBand[i] < -outerSlotDist[i])
          trackIdValues[i] = 1.0; // opening cut by the slot reads as outer
      }
      fnField.boolean_subtract(ringBand, outerSlotDist, pathValues, false);
    } else {
      pathValues = ringBand;
    }
  } else {
    // ---- Single-perimeter path (unchanged) ----
    zScalarArray outerWallField(n, 0.0);
    for (int i = 0; i < n; i++) {
      const bool isInside = (polyField[i] <= 0.0);
      double maxOuter = -(std::numeric_limits<double>::max)();
      for (int k : outerIndices) {
        const zScalar w = segmentGraphs[k].width;
        const zScalar d = segDistances[k][i];
        const zScalar phi1 = (isInside ? -d : d) + 0.5 * w;
        if (phi1 > maxOuter)
          maxOuter = phi1;
      }
      outerWallField[i] = maxOuter;
    }

    zScalarArray outerSlotDist(n, (std::numeric_limits<zScalar>::max)());
    bool hasOuterSlot = false;

    for (int k : outerIndices) {
      const zScalar trimRatio = segmentGraphs[k].trim;
      if (trimRatio > 1e-4 && trimRatio < 1.0 - 1e-4) {
        zObjGraph slotGraph;
        const zScalar w = segmentGraphs[k].width;
        if (_createPerpendicularSlot(segmentGraphs[k].graph, slotGraph,
                                     trimRatio, w * 1.5)) {
          allSlotGraphs.push_back(slotGraph);
          zScalarArray tempSlotDist;
          const float slotHalfWidth = static_cast<float>(0.5 * w);
          _evaluateGraphDistance(points, n, slotGraph, slotHalfWidth,
                                 fastSlotEdges, tempSlotDist);

          if (!hasOuterSlot) {
            outerSlotDist = tempSlotDist;
            hasOuterSlot = true;
          } else {
            for (int i = 0; i < n; i++) {
              if (tempSlotDist[i] < outerSlotDist[i])
                outerSlotDist[i] = tempSlotDist[i];
            }
          }
        }
      }
    }

    zScalarArray outerField = outerWallField;
    if (hasOuterSlot) {
      fnField.boolean_subtract(outerWallField, outerSlotDist, outerField,
                               false);
    }

    zScalarArray innerRibsField(n, (std::numeric_limits<zScalar>::max)());
    bool hasValidRibs = false;

    for (int k : ribIndices) {
      const zScalar w = segmentGraphs[k].width;
      const zScalar ribBeadW = 0.5 * w;
      zScalarArray singleRibField(n, 0.0);
      for (int i = 0; i < n; i++) {
        singleRibField[i] = segDistances[k][i] - 0.5 * ribBeadW;
      }

      const zScalar trimRatio = segmentGraphs[k].trim;
      zScalarArray singleRibTrimmed = singleRibField;
      if (trimRatio > 1e-4 && trimRatio < 1.0 - 1e-4) {
        zObjGraph slotGraph;
        if (_createPerpendicularSlot(segmentGraphs[k].graph, slotGraph,
                                     trimRatio, w * 1.5)) {
          allSlotGraphs.push_back(slotGraph);
          zScalarArray tempSlotDist;
          const float slotHalfWidth = static_cast<float>(0.5 * ribBeadW);
          _evaluateGraphDistance(points, n, slotGraph, slotHalfWidth,
                                 fastSlotEdges, tempSlotDist);
          fnField.boolean_subtract(singleRibField, tempSlotDist,
                                   singleRibTrimmed, false);
        }
      }

      if (!hasValidRibs) {
        innerRibsField = singleRibTrimmed;
        hasValidRibs = true;
      } else {
        for (int i = 0; i < n; i++) {
          if (singleRibTrimmed[i] < innerRibsField[i]) {
            innerRibsField[i] = singleRibTrimmed[i];
          }
        }
      }
    }

    if (hasValidRibs) {
      fnField.boolean_subtract(outerField, innerRibsField, pathValues, false);
    } else {
      pathValues = outerField;
    }
    // trackIdValues stays 0.0 everywhere - no inner/outer distinction here.
  }

  if (outTrimSlotsGraph && !allSlotGraphs.empty()) {
    combineGraphObjects(allSlotGraphs, *outTrimSlotsGraph);
  }

  const zScalar defaultW =
      (!segmentGraphs.empty()) ? segmentGraphs[0].width : 0.4;
  const double sdfWidth = 0.5 * (defaultW > 0 ? defaultW : 0.4);
  fnField.setFieldValues(pathValues, zFieldSDF, sdfWidth);

  if (outPathField)
    *outPathField = pathValues;
  if (outTrackIdField)
    *outTrackIdField = trackIdValues;
}

// Extracts the zero-isocontour of pathVals as raw position/edge arrays
// (marching triangles), interpolating widths and the outer/inner track tag
// onto each new crossing point. Left unsplit/uncleaned - _extractTrackPaths
// does the spur-pruning, gap-welding and inner/outer split downstream.
void _extractRawIsoContour(zObjMeshScalarField &pathFieldMesh,
                           const zScalarArray &pathVals,
                           const zScalarArray &widthVals,
                           const zScalarArray &trackIdVals,
                           zPointArray &outPositions,
                           zIntArray &outEdgeConnects,
                           zFloatArray &outVertexWeights,
                           zFloatArray &outVertexTrackIds) {
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

// Ports the spur-pruning + gap-welding + chain/loop extraction that used to
// live in zPrintSyn::buildSequencePrintPath. Splits the raw iso-contour
// first - edges whose endpoints disagree on trackId (inner vs outer) are
// dropped before cleanup, so each side is cleaned independently - then runs
// the same cleanup on each side and merges every surviving loop/chain on
// that side into one PrintPath. Callers with no track-id distinction (empty
// rawTrackIds, or an all-inner field) get a single "inner" PrintPath back.
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
    return hasTrackIds && id < int(rawTrackIds.size()) && rawTrackIds[id] >= 0.5f;
  };

  bool hasAnyOuter = false;
  if (hasTrackIds) {
    for (float v : rawTrackIds) {
      if (v >= 0.5f) {
        hasAnyOuter = true;
        break;
      }
    }
  }

  // Cleans+extracts every loop/chain whose edges match wantOuter, merging
  // all survivors into a single PrintPath (skipped if none survive).
  auto processSide = [&](bool wantOuter) {
    std::vector<std::vector<int>> adj(numV);
    std::unordered_set<unsigned long long> edgeSet;
    for (size_t e = 0; e + 1 < rawEdges.size(); e += 2) {
      int a = rawEdges[e], b = rawEdges[e + 1];
      if (a < 0 || b < 0 || a >= numV || b >= numV || a == b)
        continue;
      if (isOuterVertex(a) != wantOuter || isOuterVertex(b) != wantOuter)
        continue; // edge straddles the inner/outer boundary - drop it
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
              branchLen += vecSub(rawPositions[next], rawPositions[curr]).length();
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

    for (int i = 0; i < numV; i++) {
      if (activeAdj[i].size() == 1) {
        int next = *activeAdj[i].begin();
        if (visitedEdges.count(makeEdgeKey(i, next)))
          continue;

        std::vector<int> chain = {i};
        int prev = i;
        int curr = next;
        visitedEdges.insert(makeEdgeKey(prev, curr));

        while (true) {
          chain.push_back(curr);
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
        openChains.push_back(chain);
      }
    }

    for (int i = 0; i < numV; i++) {
      if (activeAdj[i].size() == 2) {
        for (int nb : activeAdj[i]) {
          if (visitedEdges.count(makeEdgeKey(i, nb)))
            continue;

          std::vector<int> loop = {i};
          int prev = i;
          int curr = nb;
          visitedEdges.insert(makeEdgeKey(prev, curr));
          bool loopClosed = false;

          while (true) {
            if (curr == i) {
              loopClosed = true;
              break;
            }
            loop.push_back(curr);
            if (activeAdj[curr].size() != 2)
              break;

            int nxt = -1;
            for (int n : activeAdj[curr]) {
              if (n != prev) {
                nxt = n;
                break;
              }
            }
            if (nxt < 0 || visitedEdges.count(makeEdgeKey(curr, nxt)))
              break;
            visitedEdges.insert(makeEdgeKey(curr, nxt));
            prev = curr;
            curr = nxt;
          }

          if (loopClosed && loop.size() >= 3) {
            closedLoops.push_back(loop);
          } else if (!loop.empty()) {
            openChains.push_back(loop);
          }
        }
      }
    }

    const double weldTolSq = weldTolerance * weldTolerance;
    for (auto it = openChains.begin(); it != openChains.end();) {
      if (it->size() >= 3 &&
          vecDistSq(rawPositions[it->front()], rawPositions[it->back()]) <=
              weldTolSq) {
        closedLoops.push_back(*it);
        it = openChains.erase(it);
        continue;
      }
      ++it;
    }

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

          double d0 = vecDistSq(rawPositions[c1.back()], rawPositions[c2.front()]);
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
          double d2 = vecDistSq(rawPositions[c1.front()], rawPositions[c2.back()]);
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

        if (c1.size() >= 3 &&
            vecDistSq(rawPositions[c1.front()], rawPositions[c1.back()]) <=
                weldTolSq) {
          closedLoops.push_back(c1);
          openChains.erase(openChains.begin() + bestI);
        }
      }
    }

    if (closedLoops.empty() && openChains.empty())
      return;

    zPointArray mergedPositions;
    zIntArray mergedEdges;
    zFloatArray mergedWidths;

    auto appendSeq = [&](const std::vector<int> &seq, bool isLoop) {
      const int base = static_cast<int>(mergedPositions.size());
      for (int id : seq) {
        mergedPositions.push_back(rawPositions[id]);
        mergedWidths.push_back(id < int(rawWidths.size()) ? rawWidths[id]
                                                           : 0.04f);
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
    p.isOuterTrack = wantOuter;
    zFnGraph fnP(p.graph);
    fnP.create(mergedPositions, mergedEdges);
    p.widths = mergedWidths;
    outPaths.push_back(p);
  };

  processSide(false); // inner (or the whole path, if untracked)
  if (hasAnyOuter)
    processSide(true); // outer
}

// Computes a cell-aligned bounding rectangle and resolution for the flat
// meshes.
bool _getboundingrect(const zObjMeshArray &flatMeshes,
                      zDomain<zPoint> &outBounds, double scaleFactor = 1.00,
                      double targetCellSize = 0.05, int *outResX = nullptr,
                      int *outResY = nullptr) {
  if (flatMeshes.empty())
    return false;

  zPoint minBB((std::numeric_limits<double>::max)(),
               (std::numeric_limits<double>::max)(), 0.0);
  zPoint maxBB(-(std::numeric_limits<double>::max)(),
               -(std::numeric_limits<double>::max)(), 0.0);
  bool hasValidBounds = false;

  for (const auto &mesh : flatMeshes) {
    zFnMesh fn(const_cast<zObjMesh &>(mesh));
    if (fn.numVertices() == 0)
      continue;
    zPoint meshMin, meshMax;
    fn.getBounds(meshMin, meshMax);
    minBB.x = (std::min)(minBB.x, meshMin.x);
    minBB.y = (std::min)(minBB.y, meshMin.y);
    maxBB.x = (std::max)(maxBB.x, meshMax.x);
    maxBB.y = (std::max)(maxBB.y, meshMax.y);
    hasValidBounds = true;
  }

  if (!hasValidBounds || minBB.x >= maxBB.x || minBB.y >= maxBB.y)
    return false;

  const double centerX = (minBB.x + maxBB.x) * 0.5;
  const double centerY = (minBB.y + maxBB.y) * 0.5;
  const double rawWidth = (maxBB.x - minBB.x) * scaleFactor;
  const double rawHeight = (maxBB.y - minBB.y) * scaleFactor;

  const double cellSize = (targetCellSize > 1e-6) ? targetCellSize : 0.1;
  const int resX =
      (std::max)(1, static_cast<int>(std::ceil(rawWidth / cellSize)));
  const int resY =
      (std::max)(1, static_cast<int>(std::ceil(rawHeight / cellSize)));

  const double halfW = (resX * cellSize) * 0.5;
  const double halfH = (resY * cellSize) * 0.5;

  outBounds = zDomain<zPoint>(zPoint(centerX - halfW, centerY - halfH, 0.0),
                              zPoint(centerX + halfW, centerY + halfH, 0.0));

  if (outResX)
    *outResX = resX;
  if (outResY)
    *outResY = resY;
  return true;
}

} // anonymous namespace

//--------------------------
//---- PUBLIC METHODS: STEP 3 (2D SDF SCALAR FIELD)
//--------------------------

// Builds the SDF path/width fields per layer and extracts flat 2D print
// paths, split per layer into (up to) an inner-track and an outer-track
// PrintPath.
ZSPACE_INLINE bool
zSlicer::computeField(std::vector<std::vector<PrintPath>> &outPrintPaths2D) {
  pathFields.clear();
  widthFields.clear();
  trimSlotGraphs.clear();
  outPrintPaths2D.clear();

  const int totalLayers = int(sectionMeshes2D.size());
  if (totalLayers == 0 || sectionGraphs2D.empty()) {
    m_Status =
        "failed: No unrolled 2D meshes available. Call computeUnroll() first.";
    printError("[computeField] " + m_Status);
    return false;
  }

  pathFields.assign(totalLayers, zScalarArray());
  widthFields.assign(totalLayers, zScalarArray());
  trimSlotGraphs.assign(totalLayers, zObjGraph());
  outPrintPaths2D.assign(totalLayers, std::vector<PrintPath>());

  const float printBracingWidth = i_Settings.BracingWidth;
  const float printOverlapWidth = i_Settings.OverlapWidth;
  const float targetCellSize =
      (i_Settings.FieldCellSize > 1e-6f) ? i_Settings.FieldCellSize : 0.05f;

  zDomain<zPoint> layerFieldBB;
  int resX = 0, resY = 0;
  if (!_getboundingrect(sectionMeshes2D, layerFieldBB, i_Settings.FieldScaleUp,
                        targetCellSize, &resX, &resY)) {
    m_Status =
        "failed: Could not compute bounding rectangle for unrolled meshes.";
    printError("[computeField] " + m_Status);
    return false;
  }

  const std::vector<float> &edgeWidths = i_Settings.BoundaryWidth;
  const float defaultEdgeWidth =
      (!edgeWidths.empty()) ? edgeWidths.back() : 0.5f;

  // Initialize reusable scalar field topology
  {
    zFnMeshScalarField fnField(m_ScalarField);
    fnField.clear();
    fnField.create(layerFieldBB.min, layerFieldBB.max, resX, resY, 1, true,
                   false);
  }

  const auto wallStart = std::chrono::high_resolution_clock::now();

  // Multi-threaded layer-by-layer computation (OpenMP)
#pragma omp parallel
  {
    zObjMeshScalarField localField = m_ScalarField;

#pragma omp for schedule(dynamic)
    for (int l = 0; l < totalLayers; l++) {
      std::vector<zObjGraph> boundarySegments;
      zTsUtils::segmentGraph(sectionGraphs2D[l], m_BoundarySplitIndices,
                             boundarySegments);

      zObjGraph flatBracingGraph;
      bool hasValidBracing = false;
      if (l < int(bracingGraphs2D.size())) {
        flatBracingGraph = bracingGraphs2D[l];
        zFnGraph fnBracingCheck(flatBracingGraph);
        if (fnBracingCheck.numVertices() > 0 && fnBracingCheck.numEdges() > 0) {
          hasValidBracing = true;
        }
      }

      std::vector<SegmentGraph> segmentGraphs;
      const int numBoundarySegs = int(boundarySegments.size());
      segmentGraphs.reserve(numBoundarySegs + (hasValidBracing ? 1 : 0));

      for (int seg = 0; seg < numBoundarySegs; seg++) {
        const float rawW =
            (seg < (int)edgeWidths.size()) ? edgeWidths[seg] : defaultEdgeWidth;
        const float setW = (std::max)(1e-4f, rawW);
        float boundaryTrimRatio = 0.0f;
        if (seg == i_Settings.BoundarySlotID) {
          boundaryTrimRatio = (l % 2 == 0) ? i_Settings.BoundarySlotEven
                                           : i_Settings.BoundarySlotOdd;
        }
        segmentGraphs.emplace_back(boundarySegments[seg], 0, setW, 0.5f * setW,
                                   boundaryTrimRatio);
      }

      if (hasValidBracing) {
        const float bracingTrimRatio = (l % 2 == 0) ? i_Settings.BracingSlotEven
                                                    : i_Settings.BracingSlotOdd;
        segmentGraphs.emplace_back(flatBracingGraph, 1, printBracingWidth, 0.0f,
                                   bracingTrimRatio);
      }

      std::vector<zScalarArray> segDistances;
      zObjGraph bracingSlotsGraph;
      zScalarArray widthValues;
      zScalarArray idValues;
      zScalarArray pathValues;
      zScalarArray trackIdField;

      _computeSegmentField(localField, segmentGraphs, segDistances);

      _computeWidthField(localField, segDistances, segmentGraphs,
                         i_Settings.DoubleBoundary, printOverlapWidth,
                         &widthValues, &idValues);

      _computePathField(localField, sectionGraphs2D[l], segmentGraphs,
                        segDistances, i_Settings.DoubleBoundary,
                        printOverlapWidth, &bracingSlotsGraph, &pathValues,
                        &trackIdField);

      pathFields[l] = pathValues;
      widthFields[l] = widthValues;
      trimSlotGraphs[l] = bracingSlotsGraph;

      zPointArray rawPositions;
      zIntArray rawEdges;
      zFloatArray rawWidths;
      zFloatArray rawTrackIds;
      _extractRawIsoContour(localField, pathValues, widthValues, trackIdField,
                            rawPositions, rawEdges, rawWidths, rawTrackIds);
      _extractTrackPaths(rawPositions, rawEdges, rawWidths, rawTrackIds, 0.08,
                        1.0, outPrintPaths2D[l]);
    }
  }

  const auto wallEnd = std::chrono::high_resolution_clock::now();
  const double wallTotalMs =
      std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();

  std::ostringstream out;
  out << "Compute " << totalLayers << " layers (CellSize: " << std::fixed
      << std::setprecision(2) << targetCellSize << ", Res: " << resX << "x"
      << resY << ") in " << std::setprecision(2) << wallTotalMs << " ms";
  m_Status = out.str();
  std::cout << std::endl;
  printSuccess("[computeField] " + m_Status);
  return true;
}

} // namespace zSpace
