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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_set>
#include <utility>
#include <vector>

namespace zSpace {

using namespace zTsUtils;

namespace {

inline zColor _hslToRgb(float h, float s, float l) {
  float c = (1.0f - std::abs(2.0f * l - 1.0f)) * s;
  float hp = std::fmod(h, 360.0f) / 60.0f;
  if (hp < 0.0f)
    hp += 6.0f;
  float x = c * (1.0f - std::abs(std::fmod(hp, 2.0f) - 1.0f));
  float r1 = 0.0f, g1 = 0.0f, b1 = 0.0f;
  if (hp >= 0.0f && hp < 1.0f) {
    r1 = c;
    g1 = x;
    b1 = 0.0f;
  } else if (hp >= 1.0f && hp < 2.0f) {
    r1 = x;
    g1 = c;
    b1 = 0.0f;
  } else if (hp >= 2.0f && hp < 3.0f) {
    r1 = 0.0f;
    g1 = c;
    b1 = x;
  } else if (hp >= 3.0f && hp < 4.0f) {
    r1 = 0.0f;
    g1 = x;
    b1 = c;
  } else if (hp >= 4.0f && hp < 5.0f) {
    r1 = x;
    g1 = 0.0f;
    b1 = c;
  } else {
    r1 = c;
    g1 = 0.0f;
    b1 = x;
  }
  float m = l - c * 0.5f;
  return zColor(r1 + m, g1 + m, b1 + m, 1.0f);
}

//---- Multi-Object Wavefront OBJ Exporters for Mesh Array and Graph Array

inline bool _exportMeshesToMultiOBJ(const zObjMeshArray &meshes,
                                    const std::string &filePath) {
  if (meshes.empty())
    return false;

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  ofs << "# Exported from zTsSDFSlicer (Multi-Object Mesh OBJ)\n";
  size_t globalVertexBase = 1;

  for (size_t l = 0; l < meshes.size(); l++) {
    zFnMesh fnMesh(const_cast<zObjMesh &>(meshes[l]));
    zPointArray positions;
    fnMesh.getVertexPositions(positions);
    if (positions.empty())
      continue;

    zIntArray polyConnects, polyCounts;
    fnMesh.getPolygonData(polyConnects, polyCounts);

    ofs << "o layer_" << l << "\n";
    ofs << "g layer_" << l << "\n";

    for (const auto &p : positions) {
      ofs << "v " << p.x << " " << p.y << " " << p.z << "\n";
    }

    size_t offset = 0;
    for (int count : polyCounts) {
      ofs << "f";
      for (int i = 0; i < count; i++) {
        ofs << " " << (globalVertexBase + polyConnects[offset + i]);
      }
      ofs << "\n";
      offset += count;
    }

    globalVertexBase += positions.size();
  }

  ofs.close();
  return true;
}

inline bool
_exportGraphsToMultiOBJ(const zObjGraphArray &graphs,
                        const std::string &filePath,
                        const std::string &objectPrefix = "layer_") {
  if (graphs.empty())
    return false;

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  ofs << "# Exported from zTsSDFSlicer (Multi-Object Toolpath Graph OBJ - Continuous Polyline Loops)\n";
  size_t globalVertexBase = 1;

  for (size_t l = 0; l < graphs.size(); l++) {
    zFnGraph fnGraph(const_cast<zObjGraph &>(graphs[l]));
    zPointArray positions;
    zIntArray edgeConnects;
    fnGraph.getVertexPositions(positions);
    fnGraph.getEdgeData(edgeConnects);
    const int numV = static_cast<int>(positions.size());
    const int numE = static_cast<int>(edgeConnects.size() / 2);
    if (numV < 2 || numE < 1)
      continue;

    ofs << "o " << objectPrefix << l << "\n";
    ofs << "g " << objectPrefix << l << "\n";

    // 1. Output all vertex coordinates
    for (const auto &p : positions) {
      ofs << "v " << p.x << " " << p.y << " " << p.z << "\n";
    }

    // 2. Build adjacency list to trace continuous polyline loops
    std::vector<std::vector<int>> adj(numV);
    for (int e = 0; e < numE; e++) {
      int u = edgeConnects[e * 2], v = edgeConnects[e * 2 + 1];
      if (u >= 0 && v >= 0 && u < numV && v < numV && u != v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
      }
    }

    // Trace all continuous chains and loops
    std::vector<std::vector<int>> polylineChains;
    std::unordered_set<unsigned long long> visitedEdges;

    // 2.1 Trace open chains starting from endpoints (degree 1)
    for (int i = 0; i < numV; i++) {
      if (adj[i].size() == 1) {
        int next = adj[i][0];
        if (visitedEdges.count(makeEdgeKey(i, next)))
          continue;

        std::vector<int> chain = {i};
        int prev = i;
        int curr = next;
        visitedEdges.insert(makeEdgeKey(prev, curr));

        while (true) {
          chain.push_back(curr);
          if (adj[curr].size() != 2)
            break;

          int nxt = (adj[curr][0] != prev) ? adj[curr][0] : adj[curr][1];
          if (visitedEdges.count(makeEdgeKey(curr, nxt)))
            break;
          visitedEdges.insert(makeEdgeKey(curr, nxt));
          prev = curr;
          curr = nxt;
        }
        polylineChains.push_back(chain);
      }
    }

    // 2.2 Trace closed loops (degree 2)
    for (int i = 0; i < numV; i++) {
      if (adj[i].size() == 2) {
        for (int nb : adj[i]) {
          if (visitedEdges.count(makeEdgeKey(i, nb)))
            continue;

          std::vector<int> loop = {i};
          int prev = i;
          int curr = nb;
          visitedEdges.insert(makeEdgeKey(prev, curr));
          bool isClosed = false;

          while (true) {
            if (curr == i) {
              isClosed = true;
              break;
            }
            loop.push_back(curr);
            if (adj[curr].size() != 2)
              break;

            int nxt = (adj[curr][0] != prev) ? adj[curr][0] : adj[curr][1];
            if (visitedEdges.count(makeEdgeKey(curr, nxt)))
              break;
            visitedEdges.insert(makeEdgeKey(curr, nxt));
            prev = curr;
            curr = nxt;
          }

          if (isClosed) {
            loop.push_back(i); // Explicitly close the loop
          }
          if (loop.size() >= 2) {
            polylineChains.push_back(loop);
          }
        }
      }
    }

    // 2.3 Fallback for isolated edges
    for (int e = 0; e < numE; e++) {
      int u = edgeConnects[e * 2], v = edgeConnects[e * 2 + 1];
      if (u >= 0 && v >= 0 && u < numV && v < numV &&
          !visitedEdges.count(makeEdgeKey(u, v))) {
        polylineChains.push_back({u, v});
        visitedEdges.insert(makeEdgeKey(u, v));
      }
    }

    // 3. Output each continuous polyline/loop in OBJ format as a single 'l' line
    for (const auto &chain : polylineChains) {
      if (chain.size() < 2)
        continue;
      ofs << "l";
      for (int vIdx : chain) {
        ofs << " " << (globalVertexBase + vIdx);
      }
      ofs << "\n";
    }

    globalVertexBase += positions.size();
  }

  ofs.close();
  return true;
}

inline bool _exportColoredFieldMeshesToMultiOBJ(
    const zObjMesh &baseFieldMesh,
    const std::vector<zScalarArray> &fieldScalars, const std::string &filePath,
    const std::string &objectPrefix = "field_layer_",
    bool isWidthField = false) {
  if (fieldScalars.empty())
    return false;

  zFnMesh fnBase(const_cast<zObjMesh &>(baseFieldMesh));
  zPointArray positions;
  fnBase.getVertexPositions(positions);
  if (positions.empty())
    return false;

  zIntArray polyConnects, polyCounts;
  fnBase.getPolygonData(polyConnects, polyCounts);
  if (polyCounts.empty())
    return false;

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  ofs << "# Exported Colored Scalar Fields from zTsSDFSlicer (Multi-Layer OBJ "
         "with Vertex Colors)\n";
  size_t globalVertexBase = 1;

  for (size_t l = 0; l < fieldScalars.size(); l++) {
    const auto &scalars = fieldScalars[l];
    if (scalars.empty())
      continue;

    std::vector<zColor> colors;
    if (isWidthField) {
      zTsSDFSlicer::colorWidthField(scalars, colors);
    } else {
      zTsSDFSlicer::colorPathField(scalars, colors);
    }

    ofs << "o " << objectPrefix << l << "\n";
    ofs << "g " << objectPrefix << l << "\n";

    for (size_t i = 0; i < positions.size(); i++) {
      const auto &p = positions[i];
      if (i < colors.size()) {
        const auto &c = colors[i];
        ofs << "v " << p.x << " " << p.y << " " << p.z << " " << c.r << " "
            << c.g << " " << c.b << "\n";
      } else {
        ofs << "v " << p.x << " " << p.y << " " << p.z << "\n";
      }
    }

    size_t offset = 0;
    for (int count : polyCounts) {
      ofs << "f";
      for (int i = 0; i < count; i++) {
        ofs << " " << (globalVertexBase + polyConnects[offset + i]);
      }
      ofs << "\n";
      offset += count;
    }

    globalVertexBase += positions.size();
  }

  ofs.close();
  return true;
}

} // anonymous namespace

//---- SCALAR FIELD COLORING METHODS (zTsSDFSlicer)

ZSPACE_INLINE void
zTsSDFSlicer::colorPathField(const zScalarArray &values,
                             std::vector<zColor> &outColors) {
  const int n = static_cast<int>(values.size());
  outColors.resize(n);
  if (n == 0)
    return;

  double minVal = (std::numeric_limits<double>::max)();
  double maxVal = -(std::numeric_limits<double>::max)();
  bool hasNegative = false;

  for (double value : values) {
    if (std::isfinite(value)) {
      minVal = (std::min)(minVal, value);
      maxVal = (std::max)(maxVal, value);
      if (value < 0.0)
        hasNegative = true;
    }
  }

  const zColor darkBlue(0.0f, 40.0f / 255.0f, 240.0f / 255.0f, 1.0f);
  const zColor lightBlue(180.0f / 255.0f, 200.0f / 255.0f, 1.0f, 1.0f);
  const zColor gray(220.0f / 255.0f, 220.0f / 255.0f, 220.0f / 255.0f, 1.0f);
  const zColor purple(0.45f, 0.10f, 0.85f, 1.0f);
  const zColor yellow(1.00f, 0.85f, 0.15f, 1.0f);

  for (int i = 0; i < n; i++) {
    const double value = values[i];
    if (!std::isfinite(value)) {
      outColors[i] = gray;
      continue;
    }

    if (hasNegative) {
      if (value < 0.0) {
        double t = (minVal < 0.0) ? (value / minVal) : 0.0;
        t = (std::max)(0.0, (std::min)(1.0, t));
        outColors[i] =
            zColor(float(darkBlue.r + (lightBlue.r - darkBlue.r) * t),
                   float(darkBlue.g + (lightBlue.g - darkBlue.g) * t),
                   float(darkBlue.b + (lightBlue.b - darkBlue.b) * t), 1.0f);
      } else {
        outColors[i] = gray;
      }
    } else {
      const double range = maxVal - minVal;
      double t = (range > 1e-6) ? ((value - minVal) / range) : 0.5;
      t = (std::max)(0.0, (std::min)(1.0, t));
      outColors[i] = zColor(float(purple.r + (yellow.r - purple.r) * t),
                            float(purple.g + (yellow.g - purple.g) * t),
                            float(purple.b + (yellow.b - purple.b) * t), 1.0f);
    }
  }
}

ZSPACE_INLINE void zTsSDFSlicer::colorPathField(zObjMesh &mesh,
                                                const zScalarArray &values) {
  zFnMesh fnMesh(mesh);
  zColor *colors = fnMesh.getRawVertexColors();
  if (!colors)
    return;

  std::vector<zColor> computedColors;
  zTsSDFSlicer::colorPathField(values, computedColors);

  const int count =
      (std::min)(fnMesh.numVertices(), static_cast<int>(computedColors.size()));
  for (int i = 0; i < count; i++) {
    colors[i] = computedColors[i];
  }
}

ZSPACE_INLINE zColor zTsSDFSlicer::hslToRgb(float h, float s, float l) {
  return _hslToRgb(h, s, l);
}

ZSPACE_INLINE void
zTsSDFSlicer::generateColorPalette(int numColors, std::vector<zColor> &outColors,
                                   unsigned int seed, float saturation,
                                   float lightness) {
  outColors.clear();
  if (numColors <= 0)
    return;
  outColors.reserve(numColors);
  const float step = 360.0f / static_cast<float>(numColors);
  for (int k = 0; k < numColors; k++) {
    float h = static_cast<float>(k) * step;
    outColors.push_back(zTsSDFSlicer::hslToRgb(h, saturation, lightness));
  }
  if (numColors > 1) {
    std::mt19937 rng(seed);
    std::shuffle(outColors.begin(), outColors.end(), rng);
  }
}

ZSPACE_INLINE void
zTsSDFSlicer::colorWidthField(const zScalarArray &values,
                              std::vector<zColor> &outColors) {
  const int n = static_cast<int>(values.size());
  outColors.resize(n);
  if (n == 0)
    return;

  std::vector<double> uniqueWidths;
  for (double v : values) {
    if (!std::isfinite(v) || v <= 1e-6)
      continue;
    bool found = false;
    for (double u : uniqueWidths) {
      if (std::abs(u - v) < 1e-4) {
        found = true;
        break;
      }
    }
    if (!found) {
      uniqueWidths.push_back(v);
    }
  }
  std::sort(uniqueWidths.begin(), uniqueWidths.end());
  const size_t K = uniqueWidths.size();

  std::vector<zColor> dynamicPalette;
  if (K > 0) {
    zTsSDFSlicer::generateColorPalette(static_cast<int>(K), dynamicPalette,
                                       1337, 0.85f, 0.50f);
  }

  const zColor gray(220.0f / 255.0f, 220.0f / 255.0f, 220.0f / 255.0f, 1.0f);

  for (int i = 0; i < n; i++) {
    const double value = values[i];
    if (!std::isfinite(value) || value <= 1e-6 || uniqueWidths.empty()) {
      outColors[i] = gray;
      continue;
    }

    int groupIdx = 0;
    double bestDiff = (std::numeric_limits<double>::max)();
    for (size_t g = 0; g < K; g++) {
      double diff = std::abs(uniqueWidths[g] - value);
      if (diff < bestDiff) {
        bestDiff = diff;
        groupIdx = static_cast<int>(g);
      }
    }

    outColors[i] = dynamicPalette[groupIdx % dynamicPalette.size()];
  }
}

ZSPACE_INLINE void zTsSDFSlicer::colorWidthField(zObjMesh &mesh,
                                                 const zScalarArray &values) {
  zFnMesh fnMesh(mesh);
  zColor *colors = fnMesh.getRawVertexColors();
  if (!colors)
    return;

  std::vector<zColor> computedColors;
  zTsSDFSlicer::colorWidthField(values, computedColors);

  const int count =
      (std::min)(fnMesh.numVertices(), static_cast<int>(computedColors.size()));
  for (int i = 0; i < count; i++) {
    colors[i] = computedColors[i];
  }
}

ZSPACE_INLINE void zTsSDFSlicer::translateMesh(zObjMesh &mesh,
                                               const zVector &offset) {
  zFnMesh fn(mesh);
  if (fn.numVertices() <= 0)
    return;
  zPointArray pts;
  fn.getVertexPositions(pts);
  for (auto &p : pts) {
    p.x += offset.x;
    p.y += offset.y;
    p.z += offset.z;
  }
  fn.setVertexPositions(pts);
}

ZSPACE_INLINE void zTsSDFSlicer::translateGraph(zObjGraph &graph,
                                                const zVector &offset) {
  zFnGraph fn(graph);
  if (fn.numVertices() <= 0)
    return;
  zPointArray pts;
  fn.getVertexPositions(pts);
  for (auto &p : pts) {
    p.x += offset.x;
    p.y += offset.y;
    p.z += offset.z;
  }
  fn.setVertexPositions(pts);
}

ZSPACE_INLINE bool
zTsSDFSlicer::segmentGraph(const zObjGraph &sourceBoundaryGraph,
                           const std::vector<int> &splitIndices,
                           std::vector<zObjGraph> &outBoundarySegments) {
  outBoundarySegments.clear();
  zFnGraph fnSource(const_cast<zObjGraph &>(sourceBoundaryGraph));
  zPointArray positions;
  fnSource.getVertexPositions(positions);

  const int N = int(positions.size());
  if (N < 2)
    return false;

  const int numSegs = int(splitIndices.size());
  if (numSegs == 0) {
    outBoundarySegments.push_back(sourceBoundaryGraph);
    return true;
  }

  outBoundarySegments.reserve(numSegs);

  for (int seg = 0; seg < numSegs; seg++) {
    const int startIdx = splitIndices[seg];
    const int endIdx =
        (seg + 1 < numSegs) ? splitIndices[seg + 1] : splitIndices[0];

    zPointArray outPositions;
    zIntArray outEdges;

    if (startIdx <= endIdx) {
      for (int k = startIdx; k <= endIdx; k++) {
        outPositions.push_back(positions[k]);
      }
    } else {
      for (int k = startIdx; k < N; k++) {
        outPositions.push_back(positions[k]);
      }
      for (int k = 0; k <= endIdx; k++) {
        outPositions.push_back(positions[k]);
      }
    }

    if (outPositions.size() >= 2) {
      for (size_t e = 0; e + 1 < outPositions.size(); e++) {
        outEdges.push_back(int(e));
        outEdges.push_back(int(e + 1));
      }

      zObjGraph segGraph;
      zFnGraph fnSeg(segGraph);
      fnSeg.create(outPositions, outEdges);
      outBoundarySegments.push_back(segGraph);
    }
  }

  return !outBoundarySegments.empty();
}

ZSPACE_INLINE void
zTsSDFSlicer::drawPathMesh(int layerIdx, zObjMesh &outMesh) const {
  zFnMesh fnOut(outMesh);
  fnOut.clear();

  if (layerIdx < 0 || layerIdx >= int(contourGraphs.size()))
    return;

  zFnGraph fnContour(const_cast<zObjGraph &>(contourGraphs[layerIdx]));
  zPointArray points;
  zIntArray edges;
  fnContour.getVertexPositions(points);
  fnContour.getEdgeData(edges);
  const int N = int(points.size());
  if (N < 2)
    return;

  bool closed = false;
  for (size_t e = 0; e + 1 < edges.size(); e += 2) {
    if ((edges[e] == 0 && edges[e + 1] == N - 1) ||
        (edges[e + 1] == 0 && edges[e] == N - 1)) {
      closed = true;
      break;
    }
  }

  const zFloatArray *layerWeights =
      (layerIdx < int(contourWeights.size())) ? &contourWeights[layerIdx]
                                             : nullptr;
  const zFloatArray *layerHeights =
      (layerIdx < int(contourHeights.size())) ? &contourHeights[layerIdx]
                                             : nullptr;
  const zVectorArray *layerNormals =
      (layerIdx < int(contourNormals.size())) ? &contourNormals[layerIdx]
                                             : nullptr;

  zPointArray meshPositions;
  zIntArray polyCounts;
  zIntArray polyConnects;
  meshPositions.reserve(N * 4);

  for (int i = 0; i < N; i++) {
    zPoint p = points[i];

    // 1. Tangent direction
    zVector pathDir;
    if (closed) {
      int prevIdx = (i - 1 + N) % N;
      int nextIdx = (i + 1) % N;
      pathDir = zVector(points[nextIdx].x - points[prevIdx].x,
                        points[nextIdx].y - points[prevIdx].y,
                        points[nextIdx].z - points[prevIdx].z);
    } else {
      if (i == 0) {
        pathDir = zVector(points[1].x - points[0].x,
                          points[1].y - points[0].y,
                          points[1].z - points[0].z);
      } else if (i == N - 1) {
        pathDir = zVector(points[N - 1].x - points[N - 2].x,
                          points[N - 1].y - points[N - 2].y,
                          points[N - 1].z - points[N - 2].z);
      } else {
        pathDir = zVector(points[i + 1].x - points[i - 1].x,
                          points[i + 1].y - points[i - 1].y,
                          points[i + 1].z - points[i - 1].z);
      }
    }
    if (pathDir.length() < 1e-6)
      pathDir = zVector(1, 0, 0);
    pathDir.normalize();

    // 2. 3D surface normal
    zVector norm = (layerNormals && i < int(layerNormals->size()))
                       ? (*layerNormals)[i]
                       : zVector(0, 0, 1);
    if (norm.length() < 1e-6)
      norm = zVector(0, 0, 1);
    norm.normalize();

    // 3. In-plane binormal perpendicular to tangent
    zVector binormal = norm ^ pathDir;
    if (binormal.length() < 1e-6) {
      zVector temp(0, 0, 1);
      if (std::abs(pathDir.z) > 0.9)
        temp = zVector(0, 1, 0);
      binormal = temp ^ pathDir;
    }
    binormal.normalize();

    // 4. Per-vertex width and height
    float w = (layerWeights && i < int(layerWeights->size()))
                  ? (*layerWeights)[i]
                  : 0.20f;
    if (w <= 1e-4f)
      w = 0.20f;
    float h = (layerHeights && i < int(layerHeights->size()))
                  ? (*layerHeights)[i]
                  : i_Settings.layerSpacing;
    if (h <= 1e-4f)
      h = i_Settings.layerSpacing;

    // 5. Generate 4 corner vertices for rectangular cross section
    const float halfW = w * 0.5f;
    zPoint c0(p.x + binormal.x * halfW, p.y + binormal.y * halfW,
              p.z + binormal.z * halfW); // Top-Right
    zPoint c1(p.x - binormal.x * halfW, p.y - binormal.y * halfW,
              p.z - binormal.z * halfW); // Top-Left
    zPoint c2(c1.x - norm.x * h, c1.y - norm.y * h,
              c1.z - norm.z * h); // Bottom-Left
    zPoint c3(c0.x - norm.x * h, c0.y - norm.y * h,
              c0.z - norm.z * h); // Bottom-Right

    meshPositions.push_back(c0);
    meshPositions.push_back(c1);
    meshPositions.push_back(c2);
    meshPositions.push_back(c3);
  }

  // 6. Construct 4 quads per segment along the toolpath tube
  const int segCount = closed ? N : (N - 1);
  for (int i = 0; i < segCount; i++) {
    int next = (i + 1) % N;
    // Top face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 0);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(i * 4 + 1);

    // Left face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 1);
    polyConnects.push_back(next * 4 + 1);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(i * 4 + 2);

    // Bottom face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 2);
    polyConnects.push_back(next * 4 + 2);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(i * 4 + 3);

    // Right face
    polyCounts.push_back(4);
    polyConnects.push_back(i * 4 + 3);
    polyConnects.push_back(next * 4 + 3);
    polyConnects.push_back(next * 4 + 0);
    polyConnects.push_back(i * 4 + 0);
  }

  // Cap ends for open curves
  if (!closed && N >= 2) {
    polyCounts.push_back(4);
    polyConnects.push_back(0);
    polyConnects.push_back(3);
    polyConnects.push_back(2);
    polyConnects.push_back(1);

    int last = (N - 1) * 4;
    polyCounts.push_back(4);
    polyConnects.push_back(last + 0);
    polyConnects.push_back(last + 1);
    polyConnects.push_back(last + 2);
    polyConnects.push_back(last + 3);
  }

  fnOut.create(meshPositions, polyCounts, polyConnects);
}

//---- EXPORT METHOD

ZSPACE_INLINE bool zTsSDFSlicer::exportToOBJ(const std::string &directoryPath) {
  std::string dir = directoryPath;
  if (dir.empty()) {
    const char *userProfile = std::getenv("USERPROFILE");
    dir = userProfile
              ? (std::string(userProfile) + "/Desktop/SDF_Slicer_Export")
              : "C:/Users/Default/Desktop/SDF_Slicer_Export";
  }

  std::error_code ec;
  std::filesystem::create_directories(dir, ec);

  std::vector<std::string> exportedFiles;

  // 1. 3D Section Meshes
  if (!sectionMeshes.empty()) {
    std::string path = dir + "/section_meshes.obj";
    if (_exportMeshesToMultiOBJ(sectionMeshes, path))
      exportedFiles.push_back("section_meshes.obj");
  }

  // 2. 3D Toolpath Graphs
  if (!contourGraphs.empty()) {
    std::string path = dir + "/toolpath_graphs.obj";
    if (_exportGraphsToMultiOBJ(contourGraphs, path, "toolpath_layer_"))
      exportedFiles.push_back("toolpath_graphs.obj");
  }

  // 3. Longitudinal V-Loops
  if (!m_VLoops.empty()) {
    std::string path = dir + "/vloop_graphs.obj";
    if (_exportGraphsToMultiOBJ(m_VLoops, path, "vloop_"))
      exportedFiles.push_back("vloop_graphs.obj");
  }

  // 4. 3D Bracing Graphs
  if (!bracingGraphs.empty()) {
    std::string path = dir + "/bracing_graphs.obj";
    if (_exportGraphsToMultiOBJ(bracingGraphs, path, "bracing_layer_"))
      exportedFiles.push_back("bracing_graphs.obj");
  }

  // 5. 2D Colored Width Fields (Discrete grouped colors)
  if (!widthFields.empty() && zFnMesh(m_ScalarField).numVertices() > 0) {
    std::string path = dir + "/width_fields.obj";
    if (_exportColoredFieldMeshesToMultiOBJ(m_ScalarField, widthFields, path,
                                            "width_field_layer_", true))
      exportedFiles.push_back("width_fields.obj");
  }

  // 6. 2D Colored Path Fields
  if (!pathFields.empty() && zFnMesh(m_ScalarField).numVertices() > 0) {
    std::string path = dir + "/path_fields.obj";
    if (_exportColoredFieldMeshesToMultiOBJ(m_ScalarField, pathFields, path,
                                            "path_field_layer_", false))
      exportedFiles.push_back("path_fields.obj");
  }

  // 7. 2D Toolpaths (If generated)
  if (!_2DContourGraphs.empty()) {
    std::string path = dir + "/2D_contour_graphs.obj";
    if (_exportGraphsToMultiOBJ(_2DContourGraphs, path, "2D_contour_layer_"))
      exportedFiles.push_back("2D_contour_graphs.obj");
  }

  if (exportedFiles.empty()) {
    m_Status = "failed: No geometry available to export.";
    printError("[exportToOBJ] " + m_Status);
    return false;
  }

  std::ostringstream out;
  out << "Exported " << exportedFiles.size() << " files to " << dir << ": [";
  for (size_t i = 0; i < exportedFiles.size(); i++) {
    if (i > 0)
      out << ", ";
    out << exportedFiles[i];
  }
  out << "]";
  m_Status = out.str();
  printSuccess("[exportToOBJ] " + m_Status);
  return true;
}

} // namespace zSpace
