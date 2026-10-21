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

// Parse USD/Windows SDK headers before legacy zSpace global namespace imports.
#if defined(ZSPACE_TS3DP_OPENUSD)
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/sdf/layer.h>
#include <pxr/usd/usdGeom/tokens.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/basisCurves.h>
#include <pxr/usd/usdGeom/points.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdGeom/metrics.h>
#include <pxr/base/vt/dictionary.h>
#include <pxr/base/vt/array.h>
#include <nlohmann/json.hpp>
#include <cmath>
#include <stdexcept>
#endif

#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

#ifndef ZSPACE_INLINE
#define ZSPACE_INLINE ZSPACE_TOOLSETS_INLINE
#endif

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>



namespace zSpace {
namespace zTsIO {

#if defined(ZSPACE_TS3DP_OPENUSD)
namespace {
// Only reflow the numeric arrays authored by exportPrintPathUSD. Other USD
// syntax, including strings, comments and attribute metadata, stays untouched.
std::string formatPrintPathUSDA(const std::string &text) {
  constexpr size_t maxLineLength = 500;
  std::istringstream input(text);
  std::ostringstream output;
  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t");
    const std::string declaration = first == std::string::npos ? "" : line.substr(first);
    const bool vectors = declaration.find("point3f[] points = [") == 0 ||
                         declaration.find("normal3f[] primvars:printNormals = [") == 0;
    const bool scalars = declaration.find("float[] primvars:printWidths = [") == 0 ||
                         declaration.find("float[] primvars:printHeight = [") == 0;
    const auto begin = line.find("= [");
    const auto end = line.rfind(']');
    if ((!vectors && !scalars) || begin == std::string::npos ||
        end == std::string::npos || end <= begin + 3) {
      output << line << '\n';
      continue;
    }
    const std::string indent(first, ' ');
    const std::string itemIndent = indent + "    ";
    output << line.substr(0, begin + 3) << '\n' << itemIndent;
    size_t column = itemIndent.size();
    size_t itemStart = begin + 3;
    int depth = 0;
    for (size_t i = itemStart; i <= end; ++i) {
      const char c = line[i];
      if (c == '(') ++depth;
      if (c == ')') --depth;
      if (i != end && (c != ',' || depth != 0)) continue;
      std::string item = line.substr(itemStart, i - itemStart);
      itemStart = i + 1;
      const auto start = item.find_first_not_of(" \t\r\n");
      if (start == std::string::npos) continue;
      item = item.substr(start, item.find_last_not_of(" \t\r\n") - start + 1);
      if (i != end) item += ',';
      if (column > itemIndent.size()) {
        if (column + 1 + item.size() > maxLineLength) {
          output << '\n' << itemIndent;
          column = itemIndent.size();
        } else {
          output << ' ';
          ++column;
        }
      }
      output << item;
      column += item.size();
    }
    output << '\n' << indent << line.substr(end) << '\n';
  }
  return output.str();
}
}
#endif

ZSPACE_INLINE bool exportPrintPathUSD(const PrintPath &path,
    const std::string &inputFilePath, const std::string &filePath, std::string *outError) {
  if (outError) outError->clear();
  const auto fail = [&](const std::string &message) {
    if (outError) *outError = message;
    return false;
  };
#if !defined(ZSPACE_TS3DP_OPENUSD)
  return fail("USD export requires ZSPACE_TS3DP_OPENUSD and linked OpenUSD libraries.");
#else
  try {
    const auto dot = filePath.find_last_of('.');
    std::string extension = dot == std::string::npos ? "" : filePath.substr(dot);
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return char(std::tolower(c)); });
    if (extension != ".usd" && extension != ".usda" && extension != ".usdc")
      return fail("Use a .usd, .usda or .usdc output file.");
    std::vector<zTsUtils::GraphPolyline> curves;
    if (!zTsUtils::graphPolylines(path.graph, curves) || curves.size() != 1)
      return fail("Print path must be one connected, unbranched curve.");
    zPointArray positions;
    zFnGraph(const_cast<zObjGraph &>(path.graph)).getVertexPositions(positions);
    if (path.normals.size() != positions.size() || path.widths.size() != positions.size() ||
        path.heights.size() != positions.size())
      return fail("Normals, widths and heights must have one value per graph vertex.");
    auto ids = curves.front().Vertices;
    if (ids.size() > 1 && ids.front() == ids.back()) ids.pop_back();
    if (ids.size() < (curves.front().Closed ? 3u : 2u) || ids.size() != positions.size())
      return fail("Print path contains insufficient or unused vertices.");
    pxr::VtArray<pxr::GfVec3f> points, normals;
    pxr::VtArray<float> widths, heights;
    for (int id : ids) {
      if (id < 0 || size_t(id) >= positions.size()) return fail("Invalid path vertex index.");
      const auto &p = positions[id];
      const auto &n = path.normals[id];
      if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
          !std::isfinite(n.x) || !std::isfinite(n.y) || !std::isfinite(n.z) ||
          !std::isfinite(path.widths[id]) || !std::isfinite(path.heights[id]) ||
          path.widths[id] < 0 || path.heights[id] < 0)
        return fail("Path attributes must be finite; widths and heights must be nonnegative.");
      points.push_back(pxr::GfVec3f(p.x, p.y, p.z));
      normals.push_back(pxr::GfVec3f(n.x, n.y, n.z));
      widths.push_back(path.widths[id]);
      heights.push_back(path.heights[id]);
    }
    // Load a fresh private layer so exporting over the input does not mutate
    // another live stage or reuse stale cached file contents.
    auto layer = pxr::SdfLayer::OpenAsAnonymous(inputFilePath);
    if (!layer) return fail("Could not open input USD file.");
    auto stage = pxr::UsdStage::Open(layer);
    if (!stage) return fail("Could not open input USD stage.");
    if (!stage->GetPrimAtPath(pxr::SdfPath("/block")) ||
        !stage->GetPrimAtPath(pxr::SdfPath("/bracing")) ||
        !stage->GetPrimAtPath(pxr::SdfPath("/slotpoint")))
      return fail("Input USD must contain /block, /bracing and /slotpoint.");
    if (stage->GetPrimAtPath(pxr::SdfPath("/printpath"))) {
      if (!stage->RemovePrim(pxr::SdfPath("/printpath")) ||
          stage->GetPrimAtPath(pxr::SdfPath("/printpath")))
        return fail("Cannot replace a printpath composed from another USD layer.");
    }
    auto curve = pxr::UsdGeomBasisCurves::Define(stage, pxr::SdfPath("/printpath"));
    pxr::VtArray<int> counts;
    counts.push_back(int(ids.size()));
    const bool ok = curve.CreateTypeAttr().Set(pxr::UsdGeomTokens->linear) &&
        curve.CreateWrapAttr().Set(curves.front().Closed ? pxr::UsdGeomTokens->periodic : pxr::UsdGeomTokens->nonperiodic) &&
        curve.CreateCurveVertexCountsAttr().Set(counts) &&
        curve.CreatePointsAttr().Set(points) &&
        pxr::UsdGeomPrimvarsAPI(curve.GetPrim()).CreatePrimvar(
            pxr::TfToken("printNormals"), pxr::SdfValueTypeNames->Normal3fArray,
            pxr::UsdGeomTokens->vertex).Set(normals) &&
        pxr::UsdGeomPrimvarsAPI(curve.GetPrim()).CreatePrimvar(
            pxr::TfToken("printWidths"), pxr::SdfValueTypeNames->FloatArray,
            pxr::UsdGeomTokens->vertex).Set(widths) &&
        pxr::UsdGeomPrimvarsAPI(curve.GetPrim()).CreatePrimvar(
            pxr::TfToken("printHeight"), pxr::SdfValueTypeNames->FloatArray,
            pxr::UsdGeomTokens->vertex).Set(heights);
    if (!ok) return fail("Could not author USD curve attributes.");
    if (extension == ".usda") {
      std::string text;
      if (!stage->GetRootLayer()->ExportToString(&text))
        return fail("Could not serialize USDA output.");
      const auto formatted = formatPrintPathUSDA(text);
      std::ofstream output(filePath, std::ios::binary | std::ios::trunc);
      if (!output) return fail("Could not open USDA output file.");
      output.write(formatted.data(), static_cast<std::streamsize>(formatted.size()));
      output.close();
      if (!output) return fail("Could not write USDA output file.");
    } else if (!stage->GetRootLayer()->Export(filePath)) {
      return fail("Could not write USD output file.");
    }
    return true;
  } catch (const std::exception &e) {
    return fail(e.what());
  }
#endif
}

#if defined(ZSPACE_TS3DP_OPENUSD)
namespace {
nlohmann::json usdSettingsJSON(const pxr::VtValue &value) {
  if (value.IsHolding<pxr::VtDictionary>()) {
    auto result = nlohmann::json::object();
    for (const auto &entry : value.UncheckedGet<pxr::VtDictionary>())
      result[entry.first] = usdSettingsJSON(entry.second);
    return result;
  }
  if (value.IsHolding<bool>()) return value.UncheckedGet<bool>();
  if (value.IsHolding<int>()) return value.UncheckedGet<int>();
  if (value.IsHolding<float>()) {
    const auto v = value.UncheckedGet<float>();
    if (!std::isfinite(v)) throw std::runtime_error("Settings contain a non-finite number.");
    return v;
  }
  if (value.IsHolding<double>()) {
    const auto v = value.UncheckedGet<double>();
    if (!std::isfinite(v)) throw std::runtime_error("Settings contain a non-finite number.");
    return v;
  }
  if (value.IsHolding<std::string>()) return value.UncheckedGet<std::string>();
  auto result = nlohmann::json::array();
  if (value.IsHolding<pxr::VtArray<float>>()) {
    for (auto v : value.UncheckedGet<pxr::VtArray<float>>()) result.push_back(usdSettingsJSON(pxr::VtValue(v)));
    return result;
  }
  if (value.IsHolding<pxr::VtArray<double>>()) {
    for (auto v : value.UncheckedGet<pxr::VtArray<double>>()) result.push_back(usdSettingsJSON(pxr::VtValue(v)));
    return result;
  }
  if (value.IsHolding<pxr::VtArray<int>>()) {
    for (auto v : value.UncheckedGet<pxr::VtArray<int>>()) result.push_back(v);
    return result;
  }
  throw std::runtime_error("Unsupported stage setting type; use JSON-compatible scalar or numeric array values.");
}

pxr::UsdPrim usdGeometry(const pxr::UsdStageRefPtr &stage, const std::string &path,
                         const std::string &kind) {
  const auto root = stage->GetPrimAtPath(pxr::SdfPath(path));
  if (!root) throw std::runtime_error("Missing " + kind + " object: " + path);
  if ((kind == "block" && !root.IsA<pxr::UsdGeomMesh>()) || (kind == "bracing" && !root.IsA<pxr::UsdGeomBasisCurves>() && !root.IsA<pxr::UsdGeomMesh>()) || (kind == "slotpoint" && !root.IsA<pxr::UsdGeomPoints>() && !root.IsA<pxr::UsdGeomMesh>()))
    throw std::runtime_error("Incorrect USD geometry type at " + path);
  return root;
}

zFloatArray usdWidthTable(const nlohmann::json &settings, const std::string &key) {
  const auto &table = settings.at("slicer").at(key);
  if (!table.is_array() || table.empty())
    throw std::runtime_error(key + " must be a nonempty width array.");
  zFloatArray result;
  for (const auto &entry : table) {
    if (!entry.is_number()) throw std::runtime_error(key + " contains a nonnumeric width.");
    const float width = entry.get<float>();
    if (!std::isfinite(width) || width <= 0)
      throw std::runtime_error(key + " requires finite positive widths.");
    result.push_back(width);
  }
  return result;
}

pxr::VtArray<int> usdUniformIds(const pxr::UsdPrim &prim, const char *name, size_t count, bool optional = false) {
  const auto pv = pxr::UsdGeomPrimvarsAPI(prim).GetPrimvar(pxr::TfToken(name));
  pxr::VtArray<int> values;
  // Desktop-converted files used custom arrays before the primvar schema.
  // An authored primvar takes precedence and must satisfy its own contract.
  if (pv) {
    if (pv.GetInterpolation() != pxr::UsdGeomTokens->uniform ||
        !pv.ComputeFlattened(&values) || values.size() != count)
      throw std::runtime_error(prim.GetPath().GetString() + ": " + name + " must be a uniform int[] with one value per face/curve.");
  } else {
    const auto attribute = prim.GetAttribute(pxr::TfToken(name));
    if (!attribute && optional) { values.assign(count, 0); return values; }
    if (!attribute || !attribute.Get(&values) || values.size() != count)
      throw std::runtime_error(prim.GetPath().GetString() + ": " + name + " must be an int[] with one value per face/edge.");
  }
  return values;
}
zPointArray usdWorldPoints(const pxr::UsdGeomPointBased &geometry, pxr::UsdGeomXformCache &cache) {
  pxr::VtArray<pxr::GfVec3f> points;
  if (!geometry.GetPointsAttr().Get(&points)) throw std::runtime_error("Missing points on " + geometry.GetPath().GetString());
  const auto matrix = cache.GetLocalToWorldTransform(geometry.GetPrim());
  zPointArray result;
  result.reserve(points.size());
  for (const auto &point : points) {
    const auto p = matrix.Transform(pxr::GfVec3d(point));
    if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]))
      throw std::runtime_error("Geometry contains non-finite transformed positions.");
    result.push_back(zPoint(p[0], p[1], p[2]));
  }
  return result;
}
}
#endif

ZSPACE_INLINE bool readSlicerInputUSD(const std::string &filePath, SlicerUSDInput &out,
    std::string *outError, const std::string &blockPath,
    const std::string &bracingPath, const std::string &startPath) {
  out = SlicerUSDInput();
  if (outError) outError->clear();
#if !defined(ZSPACE_TS3DP_OPENUSD)
  if (outError) *outError = "USD input requires ZSPACE_TS3DP_OPENUSD and linked OpenUSD libraries.";
  return false;
#else
  try {
    auto stage = pxr::UsdStage::Open(filePath);
    if (!stage) throw std::runtime_error("Could not open USD stage: " + filePath);
    const auto metadata = stage->GetRootLayer()->GetCustomLayerData();
    nlohmann::json settings = nlohmann::json::object();
    for (const std::string group : {"slicer", "printSyn"}) {
      const auto it = metadata.find(group);
      if (it == metadata.end() || !it->second.IsHolding<pxr::VtDictionary>())
        throw std::runtime_error("Stage customLayerData requires dictionary " + group + ".");
      settings[group] = usdSettingsJSON(it->second);
    }
    out.SettingsJSON = settings.dump();
    out.MetersPerUnit = pxr::UsdGeomGetStageMetersPerUnit(stage);
    out.UpAxis = pxr::UsdGeomGetStageUpAxis(stage).GetString();
    if (!(out.MetersPerUnit > 0) || !std::isfinite(out.MetersPerUnit))
      throw std::runtime_error("Invalid stage metersPerUnit.");
    const auto block = usdGeometry(stage, blockPath, "block");
    const auto bracing = usdGeometry(stage, bracingPath, "bracing");
    const auto start = usdGeometry(stage, startPath, "slotpoint");
    pxr::UsdGeomXformCache cache(pxr::UsdTimeCode::Default());
    pxr::UsdGeomMesh mesh(block);
    auto positions = usdWorldPoints(mesh, cache);
    pxr::VtArray<int> counts, indices, tags, holes;
    if (!mesh.GetFaceVertexCountsAttr().Get(&counts) || !mesh.GetFaceVertexIndicesAttr().Get(&indices) || counts.empty())
      throw std::runtime_error("Block has no valid polygon topology.");
    mesh.GetHoleIndicesAttr().Get(&holes);
    if (!holes.empty()) throw std::runtime_error("Block contains USD hole faces; bake a closed polygon mesh first.");
    tags = usdUniformIds(block, "FaceTag", counts.size());
    size_t total = 0;
    for (int count : counts) {
      if (count < 3) throw std::runtime_error("Block contains a face with fewer than three vertices.");
      total += static_cast<size_t>(count);
    }
    if (total != indices.size()) throw std::runtime_error("Block face counts do not match connectivity.");
    for (int id : indices) if (id < 0 || static_cast<size_t>(id) >= positions.size())
      throw std::runtime_error("Block has an out-of-range vertex index.");
    const auto boundaryWidths = usdWidthTable(settings, "BoundaryWidth");
    const auto bracingWidths = usdWidthTable(settings, "BracingWidth");
    const auto faceWidthIds = usdUniformIds(block, "widthid", counts.size());
    for (size_t f = 0; f < tags.size(); ++f) {
      if (tags[f] != Bottom && tags[f] != Top && tags[f] != Side)
        throw std::runtime_error("FaceTag must use Bottom=1, Top=2, Side=3.");
      if (tags[f] == Side) {
        if (faceWidthIds[f] < 0 || static_cast<size_t>(faceWidthIds[f]) >= boundaryWidths.size())
          throw std::runtime_error("block face widthid is outside BoundaryWidth.");
      } else if (faceWidthIds[f] != -1)
        throw std::runtime_error("top/bottom face widthid must be -1.");
    }
    out.FaceWidthIds.assign(faceWidthIds.begin(), faceWidthIds.end());
    out.FaceTags.assign(tags.begin(), tags.end());
    zIntArray polygonCounts(counts.begin(), counts.end()), connects(indices.begin(), indices.end());
    pxr::TfToken orientation;
    mesh.GetOrientationAttr().Get(&orientation);
    const bool reverse = (orientation == pxr::UsdGeomTokens->leftHanded) !=
        (cache.GetLocalToWorldTransform(block).GetDeterminant() < 0);
    if (reverse) {
      size_t offset = 0;
      for (int count : polygonCounts) { std::reverse(connects.begin()+offset, connects.begin()+offset+count); offset += count; }
    }
    zFnMesh fn(out.Block);
    fn.create(positions, polygonCounts, connects);
    if (fn.numVertices() != static_cast<int>(positions.size()) || fn.numPolygons() != static_cast<int>(counts.size()))
      throw std::runtime_error("Mesh creation changed topology; face attributes cannot be aligned.");
    if (bracing.IsA<pxr::UsdGeomMesh>()) {
      // Legacy zSpace graph encoding: Mesh points plus custom edge pairs.
      const auto points = usdWorldPoints(pxr::UsdGeomPointBased(bracing), cache);
      pxr::VtArray<int> indices;
      if (!bracing.GetAttribute(pxr::TfToken("edgeVertexIndices")).Get(&indices) ||
          indices.empty() || indices.size() % 2)
        throw std::runtime_error("Mesh bracing requires nonempty edgeVertexIndices pairs.");
      const auto widths = usdUniformIds(bracing, "widthid", indices.size() / 2, true);
      for (size_t e = 0; e < indices.size(); e += 2) {
        const int a = indices[e], b = indices[e+1];
        if (a < 0 || b < 0 || size_t(a) >= points.size() || size_t(b) >= points.size() ||
            zTsUtils::vecDistSq(points[a], points[b]) <= 0)
          throw std::runtime_error("Bracing edge contains an invalid index or zero-length segment.");
      }
      zObjGraph graph; zPointArray ownedPoints = points;
      zIntArray edges(indices.begin(), indices.end()); zFnGraph fnGraph(graph);
      fnGraph.create(ownedPoints, edges);
      if (fnGraph.numEdges() != int(indices.size()/2))
        throw std::runtime_error("Bracing graph creation changed edge order/count; widths cannot be aligned.");
      std::vector<zTsUtils::GraphPolyline> curves;
      if (!zTsUtils::graphPolylines(graph, curves))
        throw std::runtime_error("Bracing requires unbranched chains or loops.");
      out.BracingCurves.resize(curves.size());
      for (size_t c = 0; c < curves.size(); ++c) {
        const auto& curve = curves[c]; const int widthId = widths[curve.Edges.front()];
        if (widthId < 0 || size_t(widthId) >= bracingWidths.size())
          throw std::runtime_error("Bracing widthid is outside BracingWidth.");
        for (int e : curve.Edges) if (widths[e] != widthId)
          throw std::runtime_error("Each connected Mesh bracing curve must use one widthid.");
        zPointArray cp; zIntArray ce;
        for (int v : curve.Vertices) cp.push_back(points[v]);
        for (size_t j = 0; j < curve.Edges.size(); ++j) {ce.push_back(j);ce.push_back((j+1)%cp.size());}
        zFnGraph(out.BracingCurves[c]).create(cp,ce);
        out.BracingWidthIds.push_back(widthId); out.BracingPrintWidths.push_back(bracingWidths[widthId]);
      }
    } else {
      pxr::UsdGeomBasisCurves curves(bracing);
      pxr::TfToken curveType, wrap;
      curves.GetTypeAttr().Get(&curveType);
      curves.GetWrapAttr().Get(&wrap);
      if (curveType != pxr::UsdGeomTokens->linear ||
          (wrap != pxr::UsdGeomTokens->nonperiodic && wrap != pxr::UsdGeomTokens->periodic))
        throw std::runtime_error("bracing requires linear periodic or nonperiodic BasisCurves.");
      const auto curvePoints = usdWorldPoints(curves, cache);
      pxr::VtArray<int> curveCounts;
      if (!curves.GetCurveVertexCountsAttr().Get(&curveCounts) || curveCounts.empty())
        throw std::runtime_error("bracing has no curves.");
      const bool closed = wrap == pxr::UsdGeomTokens->periodic;
      size_t pointCount = 0;
      for (int count : curveCounts) {
        if (count < (closed ? 3 : 2)) throw std::runtime_error("bracing curve has too few points.");
        pointCount += size_t(count);
      }
      if (pointCount != curvePoints.size()) throw std::runtime_error("bracing curve counts do not match points.");
      const auto widthIds = usdUniformIds(bracing, "widthid", curveCounts.size(), true);
      out.BracingWidthIds.assign(widthIds.begin(), widthIds.end());
      out.BracingCurves.resize(curveCounts.size());
      size_t offset = 0;
      for (size_t c = 0; c < curveCounts.size(); ++c) {
        const int count = curveCounts[c], id = widthIds[c];
        if (id < 0 || size_t(id) >= bracingWidths.size())
          throw std::runtime_error("bracing curve widthid is outside BracingWidth.");
        zPointArray ribPoints(curvePoints.begin() + offset, curvePoints.begin() + offset + count);
        zIntArray ribEdges;
        for (int j = 0; j < (closed ? count : count - 1); ++j) {
          const int next = (j + 1) % count;
          if (zTsUtils::vecDistSq(ribPoints[j], ribPoints[next]) <= 0)
            throw std::runtime_error("bracing contains a zero-length segment.");
          ribEdges.push_back(j); ribEdges.push_back(next);
        }
        zFnGraph(out.BracingCurves[c]).create(ribPoints, ribEdges);
        out.BracingPrintWidths.push_back(bracingWidths[id]);
        offset += count;
      }
    }
    const auto anchors = usdWorldPoints(pxr::UsdGeomPointBased(start), cache);
    if (anchors.size() != 1) throw std::runtime_error("slotpoint must contain exactly one point.");
    out.StartPoint = anchors.front();
    return true;
  } catch (const std::exception &error) {
    out = SlicerUSDInput();
    if (outError) *outError = error.what();
    return false;
  }
#endif
}

} // namespace zTsIO
} // namespace zSpace
