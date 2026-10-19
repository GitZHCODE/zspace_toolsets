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

#include <headers/zToolsets/ts3dp/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>
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
                                   std::vector<zColor> &outColors,
                                   const zScalarArray *paletteWidths) {
  const int n = static_cast<int>(values.size());
  outColors.resize(n);
  if (n == 0)
    return;

  const zScalarArray &paletteSource =
      (paletteWidths && !paletteWidths->empty()) ? *paletteWidths : values;
  std::vector<double> uniqueWidths;
  for (double v : paletteSource) {
    if (!std::isfinite(v) || v <= 1e-3)
      continue;
    bool found = false;
    for (double u : uniqueWidths) {
      if (std::abs(u - v) < 0.1) {
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
    if (!std::isfinite(value) || value <= 1e-3 || uniqueWidths.empty()) {
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

ZSPACE_INLINE void colorWidthField(zObjMesh &mesh, const zScalarArray &values,
                                   const zScalarArray *paletteWidths) {
  zFnMesh fnMesh(mesh);
  zColor *colors = fnMesh.getRawVertexColors();
  if (!colors)
    return;

  std::vector<zColor> computedColors;
  colorWidthField(values, computedColors, paletteWidths);

  const int count =
      (std::min)(fnMesh.numVertices(), static_cast<int>(computedColors.size()));
  for (int i = 0; i < count; i++) {
    colors[i] = computedColors[i];
  }
}

} // namespace zTsUtils
} // namespace zSpace
