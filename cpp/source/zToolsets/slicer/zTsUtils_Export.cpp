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

#include <headers/zToolsets/slicer/zTsUtils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <unordered_set>
#include <vector>

namespace zSpace {
namespace zTsUtils {

//---- Color Format Conversion

ZSPACE_INLINE zColor hslToRgb(float h, float s, float l) {
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

//---- Distinct Shuffled Color Palette Generator

ZSPACE_INLINE void generateColorPalette(int numColors,
                                        std::vector<zColor> &outColors,
                                        unsigned int seed, float saturation,
                                        float lightness) {
  outColors.clear();
  if (numColors <= 0)
    return;
  outColors.reserve(numColors);
  const float step = 360.0f / static_cast<float>(numColors);
  for (int k = 0; k < numColors; k++) {
    float h = static_cast<float>(k) * step;
    outColors.push_back(hslToRgb(h, saturation, lightness));
  }
  if (numColors > 1) {
    std::mt19937 rng(seed);
    std::shuffle(outColors.begin(), outColors.end(), rng);
  }
}

//---- Field-Mesh Coloring: Positive Values Solid, Negative Values a Gradient

ZSPACE_INLINE void colorPathField(const zScalarArray &values,
                                  std::vector<zColor> &outColors,
                                  const zColor &positiveColor,
                                  const zColor &negativeColorNearZero,
                                  const zColor &negativeColorAtMin) {
  const int n = static_cast<int>(values.size());
  outColors.resize(n);
  if (n == 0)
    return;

  double minVal = (std::numeric_limits<double>::max)();
  for (double value : values) {
    if (std::isfinite(value)) {
      minVal = (std::min)(minVal, value);
    }
  }

  const zColor gray(220.0f / 255.0f, 220.0f / 255.0f, 220.0f / 255.0f, 1.0f);

  for (int i = 0; i < n; i++) {
    const double value = values[i];
    if (!std::isfinite(value)) {
      outColors[i] = gray;
      continue;
    }

    if (value < 0.0 && minVal < 0.0) {
      double t = value / minVal; // 0 at zero, 1 at the most negative value
      t = (std::max)(0.0, (std::min)(1.0, t));
      outColors[i] = zColor(
          float(negativeColorNearZero.r +
                (negativeColorAtMin.r - negativeColorNearZero.r) * t),
          float(negativeColorNearZero.g +
                (negativeColorAtMin.g - negativeColorNearZero.g) * t),
          float(negativeColorNearZero.b +
                (negativeColorAtMin.b - negativeColorNearZero.b) * t),
          1.0f);
    } else {
      outColors[i] = positiveColor;
    }
  }
}

ZSPACE_INLINE void colorPathField(zObjMesh &mesh, const zScalarArray &values,
                                  const zColor &positiveColor,
                                  const zColor &negativeColorNearZero,
                                  const zColor &negativeColorAtMin) {
  zFnMesh fnMesh(mesh);
  zColor *colors = fnMesh.getRawVertexColors();
  if (!colors)
    return;

  std::vector<zColor> computedColors;
  colorPathField(values, computedColors, positiveColor, negativeColorNearZero,
                 negativeColorAtMin);

  const int count =
      (std::min)(fnMesh.numVertices(), static_cast<int>(computedColors.size()));
  for (int i = 0; i < count; i++) {
    colors[i] = computedColors[i];
  }
}

//---- Field-Mesh Coloring: Discrete Index Lookup From a Generated Color List

ZSPACE_INLINE void colorWidthField(const zScalarArray &values,
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
    generateColorPalette(static_cast<int>(K), dynamicPalette, 1337, 0.85f,
                         0.50f);
  }

  const zColor fallbackColor(1.0f, 1.0f, 1.0f, 1.0f);

  for (int i = 0; i < n; i++) {
    const double value = values[i];
    if (!std::isfinite(value) || value <= 1e-6 || uniqueWidths.empty()) {
      outColors[i] = fallbackColor;
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

ZSPACE_INLINE void colorWidthField(zObjMesh &mesh, const zScalarArray &values) {
  zFnMesh fnMesh(mesh);
  zColor *colors = fnMesh.getRawVertexColors();
  if (!colors)
    return;

  std::vector<zColor> computedColors;
  colorWidthField(values, computedColors);

  const int count =
      (std::min)(fnMesh.numVertices(), static_cast<int>(computedColors.size()));
  for (int i = 0; i < count; i++) {
    colors[i] = computedColors[i];
  }
}

//---- Generic Multi-Object Wavefront OBJ Exporters

ZSPACE_INLINE bool exportMeshToOBJ(const zObjMeshArray &meshes,
                                          const std::string &filePath) {
  if (meshes.empty())
    return false;

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  ofs << "# Exported from zTsUtils (Multi-Object Mesh OBJ)\n";
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

ZSPACE_INLINE bool
exportGraphToOBJ(const zObjGraphArray &graphs, const std::string &filePath,
                       const std::string &objectPrefix) {
  if (graphs.empty())
    return false;

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  ofs << "# Exported from zTsUtils (Multi-Object Toolpath Graph OBJ - "
         "Continuous Polyline Loops)\n";
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

ZSPACE_INLINE bool exportColoredMeshToOBJ(
    const zObjMesh &baseMesh, const std::vector<std::vector<zColor>> &layerColors,
    const std::string &filePath, const std::string &objectPrefix) {
  if (layerColors.empty())
    return false;

  zFnMesh fnBase(const_cast<zObjMesh &>(baseMesh));
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

  ofs << "# Exported from zTsUtils (Multi-Layer Colored Mesh OBJ with Vertex "
         "Colors)\n";
  size_t globalVertexBase = 1;

  for (size_t l = 0; l < layerColors.size(); l++) {
    const auto &colors = layerColors[l];
    if (colors.empty())
      continue;

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

} // namespace zTsUtils
} // namespace zSpace
