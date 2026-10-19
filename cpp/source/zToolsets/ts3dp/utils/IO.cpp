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
#include <cctype>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

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
                         declaration.find("normal3f[] normals = [") == 0;
    const bool scalars = declaration.find("float[] widths = [") == 0 ||
                         declaration.find("float[] primvars:height = [") == 0;
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
    const std::string &filePath, std::string *outError) {
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
    auto stage = pxr::UsdStage::CreateInMemory();
    if (!stage) return fail("Could not create USD stage.");
    pxr::UsdGeomSetStageMetersPerUnit(stage, 0.001);
    pxr::UsdGeomSetStageUpAxis(stage, pxr::UsdGeomTokens->z);
    auto curve = pxr::UsdGeomBasisCurves::Define(stage, pxr::SdfPath("/printpath"));
    stage->SetDefaultPrim(curve.GetPrim());
    pxr::VtArray<int> counts;
    counts.push_back(int(ids.size()));
    const bool ok = curve.CreateTypeAttr().Set(pxr::UsdGeomTokens->linear) &&
        curve.CreateWrapAttr().Set(curves.front().Closed ? pxr::UsdGeomTokens->periodic : pxr::UsdGeomTokens->nonperiodic) &&
        curve.CreateCurveVertexCountsAttr().Set(counts) &&
        curve.CreatePointsAttr().Set(points) &&
        curve.CreateNormalsAttr().Set(normals) &&
        curve.SetNormalsInterpolation(pxr::UsdGeomTokens->vertex) &&
        curve.CreateWidthsAttr().Set(widths) &&
        curve.SetWidthsInterpolation(pxr::UsdGeomTokens->vertex) &&
        pxr::UsdGeomPrimvarsAPI(curve.GetPrim()).CreatePrimvar(
            pxr::TfToken("height"), pxr::SdfValueTypeNames->FloatArray,
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
  if ((kind == "block" && !root.IsA<pxr::UsdGeomMesh>()) || (kind == "bracing" && !root.IsA<pxr::UsdGeomBasisCurves>()) || (kind == "slotpoint" && !root.IsA<pxr::UsdGeomPoints>()))
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
  if (!pv && optional) { values.assign(count, 0); return values; }
  if (!pv || pv.GetInterpolation() != pxr::UsdGeomTokens->uniform ||
      !pv.ComputeFlattened(&values) || values.size() != count)
    throw std::runtime_error(prim.GetPath().GetString() + ": " + name + " must be a uniform int[] with one value per face/curve.");
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
    const auto anchors = usdWorldPoints(pxr::UsdGeomPoints(start), cache);
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

//--------------------------
//---- PRIVATE HELPERS
//--------------------------

namespace {

// Parses the rest of a "v" line: x y z, then an optional r g b. Vertices
// without a color get black so colors stay index-aligned with positions.
void _parseVertex(std::istringstream &iss, zPointArray &positions,
                  zColorArray &colors, bool &anyColor) {
  std::vector<double> values;
  double value;
  while (iss >> value)
    values.push_back(value);
  while (values.size() < 3)
    values.push_back(0.0);
  positions.push_back(zPoint(values[0], values[1], values[2]));
  if (values.size() >= 6) {
    colors.push_back(zColor(static_cast<float>(values[3]),
                            static_cast<float>(values[4]),
                            static_cast<float>(values[5]), 1.0f));
    anyColor = true;
  } else {
    colors.push_back(zColor(0.0f, 0.0f, 0.0f, 1.0f));
  }
}

// Parses the vertex index of an "f"/"l" token ("7", "7/2", "7//3", "7/2/3",
// or negative = relative to the vertices read so far) into a 0-based index;
// -1 if malformed or out of range.
int _parseIndex(const std::string &token, int vertexCount) {
  try {
    const int raw = std::stoi(token.substr(0, token.find('/')));
    const int index = (raw < 0) ? vertexCount + raw : raw - 1;
    return (index >= 0 && index < vertexCount) ? index : -1;
  } catch (...) {
    return -1;
  }
}

// Reads every index token left on an "f"/"l" line; false on any bad index.
bool _parseIndices(std::istringstream &iss, int vertexCount,
                   zIntArray &outIndices) {
  outIndices.clear();
  std::string token;
  while (iss >> token) {
    const int index = _parseIndex(token, vertexCount);
    if (index < 0)
      return false;
    outIndices.push_back(index);
  }
  return true;
}

// OBJ indices are global across objects. Compact only after parsing indices.
void _compactOBJPositions(zPointArray &positions, zIntArray &connects) {
  zIntArray localIds(positions.size(), -1);
  for (int id : connects) localIds[id] = 0;
  zPointArray compact;
  for (size_t i = 0; i < positions.size(); ++i) {
    if (localIds[i] < 0) continue;
    localIds[i] = static_cast<int>(compact.size());
    compact.push_back(positions[i]);
  }
  for (int &id : connects) id = localIds[id];
  positions.swap(compact);
}

// Positions written with enough digits to read back bit-identical.
void _writePositions(std::ofstream &ofs, const zPointArray &positions,
                     const zColorArray &colors) {
  const auto defaultPrecision = ofs.precision();
  for (size_t i = 0; i < positions.size(); i++) {
    ofs.precision(std::numeric_limits<double>::max_digits10);
    ofs << "v " << double(positions[i].x) << " " << double(positions[i].y)
        << " " << double(positions[i].z);
    if (!colors.empty()) {
      ofs.precision(defaultPrecision);
      ofs << " " << colors[i].r << " " << colors[i].g << " " << colors[i].b;
    }
    ofs << "\n";
  }
  ofs.precision(defaultPrecision);
}

} // anonymous namespace

//--------------------------
//---- READERS
//--------------------------

ZSPACE_INLINE bool parseFaceMaterial(const std::string &name, int &outTag,
                                     int &outWidthId) {
  std::string key = name;
  std::transform(key.begin(), key.end(), key.begin(),
                 [](unsigned char c) { return char(std::tolower(c)); });
  // Blender duplicates: "side_1.001" -> "side_1".
  const size_t dot = key.rfind('.');
  if (dot != std::string::npos && dot + 1 < key.size() &&
      std::all_of(key.begin() + dot + 1, key.end(),
                  [](unsigned char c) { return std::isdigit(c) != 0; }))
    key.erase(dot);

  outWidthId = 0;
  if (key == "bottom") {
    outTag = static_cast<int>(VertexTag::Bottom);
    return true;
  }
  if (key == "top") {
    outTag = static_cast<int>(VertexTag::Top);
    return true;
  }
  if (key == "side") {
    outTag = static_cast<int>(VertexTag::Side);
    return true;
  }
  if (key.rfind("side_", 0) == 0 && key.size() > 5 &&
      std::all_of(key.begin() + 5, key.end(),
                  [](unsigned char c) { return std::isdigit(c) != 0; })) {
    try {
      outWidthId = std::stoi(key.substr(5));
    } catch (...) {
      return false;
    }
    outTag = static_cast<int>(VertexTag::Side);
    return true;
  }
  return false;
}

ZSPACE_INLINE bool readMeshFromOBJ(const std::string &filePath,
                                   zObjMesh &outMesh, zIntArray &outFaceTags,
                                   zIntArray &outFaceWidthIds,
                                   std::string *outError,
                                   const std::string &objectName) {
  outFaceTags.clear();
  outFaceWidthIds.clear();
  const auto fail = [&](const std::string &reason) {
    if (outError)
      *outError = reason;
    return false;
  };
  std::ifstream ifs(filePath);
  if (!ifs.is_open())
    return fail("cannot open " + filePath);

  zPointArray positions;
  zColorArray colors; // parsed only to skip them
  bool anyColor = false;
  zIntArray polyCounts, polyConnects, face, faceTags, faceWidthIds;
  // Material of the current "usemtl" block: -1 = none yet / unknown name.
  int currentTag = -1, currentWidthId = 0;
  std::string currentMaterial;
  int untagged = 0, firstUntaggedLine = -1, lineNumber = 0;
  std::string currentObject;
  std::string line;
  while (std::getline(ifs, line)) {
    lineNumber++;
    std::istringstream iss(line);
    std::string tag;
    iss >> tag;
    if (tag == "o") {
      iss >> currentObject;
    } else if (tag == "v") {
      _parseVertex(iss, positions, colors, anyColor);
    } else if (tag == "usemtl" && (objectName.empty() || currentObject == objectName)) {
      std::getline(iss >> std::ws, currentMaterial);
      while (!currentMaterial.empty() &&
             std::isspace(static_cast<unsigned char>(currentMaterial.back())))
        currentMaterial.pop_back();
      if (!parseFaceMaterial(currentMaterial, currentTag, currentWidthId))
        return fail("line " + std::to_string(lineNumber) + ": unknown material \"" +
                    currentMaterial +
                    "\" (expected bottom, top, side or side_<n>)");
    } else if (tag == "f" && (objectName.empty() || currentObject == objectName)) {
      if (!_parseIndices(iss, int(positions.size()), face))
        return fail("line " + std::to_string(lineNumber) +
                    ": invalid face index");
      if (face.size() >= 3) {
        if (currentTag < 0) {
          if (firstUntaggedLine < 0)
            firstUntaggedLine = lineNumber;
          untagged++;
        }
        polyCounts.push_back(int(face.size()));
        polyConnects.insert(polyConnects.end(), face.begin(), face.end());
        faceTags.push_back(currentTag);
        faceWidthIds.push_back(currentWidthId);
      }
    }
  }

  if (positions.empty() || polyCounts.empty())
    return fail("no \"v\"/\"f\" lines");
  if (untagged > 0)
    return fail(std::to_string(untagged) +
                " faces have no \"usemtl\" before them (first on line " +
                std::to_string(firstUntaggedLine) + ")");

  zFnMesh fnMesh(outMesh);
  fnMesh.clear();
  if (!objectName.empty()) _compactOBJPositions(positions, polyConnects);
  fnMesh.create(positions, polyCounts, polyConnects);
  if (fnMesh.numPolygons() != int(faceTags.size()))
    return fail("mesh has " + std::to_string(fnMesh.numPolygons()) +
                " faces after creation, " + std::to_string(faceTags.size()) +
                " read - face tags would not line up");
  outFaceTags = faceTags;
  outFaceWidthIds = faceWidthIds;
  return true;
}

ZSPACE_INLINE bool readGraphFromOBJ(const std::string &filePath,
                                    zObjGraph &outGraph,
                                    const std::string &objectName) {
  std::ifstream ifs(filePath);
  if (!ifs.is_open())
    return false;

  zPointArray positions;
  zColorArray colors; // parsed only to skip them - graphs carry no colors
  bool anyColor = false;
  zIntArray edgeConnects, indices;
  std::string currentObject;
  std::string line;
  while (std::getline(ifs, line)) {
    std::istringstream iss(line);
    std::string tag;
    iss >> tag;
    if (tag == "o") {
      iss >> currentObject;
    } else if (tag == "v") {
      _parseVertex(iss, positions, colors, anyColor);
    } else if (tag == "l" && (objectName.empty() || currentObject == objectName)) {
      if (!_parseIndices(iss, int(positions.size()), indices))
        return false;
      for (size_t k = 0; k + 1 < indices.size(); k++) {
        edgeConnects.push_back(indices[k]);
        edgeConnects.push_back(indices[k + 1]);
      }
    }
  }

  if (positions.empty() || edgeConnects.empty())
    return false;

  zFnGraph fnGraph(outGraph);
  fnGraph.clear();
  if (!objectName.empty()) _compactOBJPositions(positions, edgeConnects);
  fnGraph.create(positions, edgeConnects);
  return true;
}

ZSPACE_INLINE bool readPointsFromOBJ(const std::string &filePath,
                                     zPointArray &outPoints,
                                     const std::string &objectName) {
  std::ifstream ifs(filePath);
  if (!ifs.is_open())
    return false;

  zPointArray positions;
  zColorArray colors; // parsed only to skip them
  bool anyColor = false;
  std::string currentObject;
  std::string line;
  while (std::getline(ifs, line)) {
    std::istringstream iss(line);
    std::string tag;
    iss >> tag;
    if (tag == "o")
      iss >> currentObject;
    else if (tag == "v" && (objectName.empty() || currentObject == objectName))
      _parseVertex(iss, positions, colors, anyColor);
  }

  if (positions.empty())
    return false;
  outPoints = positions;
  return true;
}

//--------------------------
//---- WRITERS
//--------------------------

ZSPACE_INLINE bool exportMeshToOBJ(const zObjMesh &mesh,
                                   const std::string &filePath,
                                   const zIntArray &faceTags,
                                   const zIntArray &faceWidthIds) {
  zFnMesh fnMesh(const_cast<zObjMesh &>(mesh));
  zPointArray positions;
  fnMesh.getVertexPositions(positions);
  zIntArray polyConnects, polyCounts;
  fnMesh.getPolygonData(polyConnects, polyCounts);
  if (!faceTags.empty() && faceTags.size() != polyCounts.size())
    return false;
  if (!faceWidthIds.empty() && faceWidthIds.size() != polyCounts.size())
    return false;

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  ofs << "# zTsIO mesh: " << positions.size() << " vertices, "
      << polyCounts.size() << " faces"
      << (faceTags.empty() ? "" : ", face tags as usemtl") << "\n";
  _writePositions(ofs, positions, zColorArray());

  const auto materialOf = [&](size_t f) {
    const int tag = faceTags[f];
    if (tag == static_cast<int>(VertexTag::Bottom))
      return std::string("bottom");
    if (tag == static_cast<int>(VertexTag::Top))
      return std::string("top");
    return "side_" + std::to_string(faceWidthIds.empty() ? 0 : faceWidthIds[f]);
  };
  std::string lastMaterial;
  size_t offset = 0;
  for (size_t f = 0; f < polyCounts.size(); f++) {
    const int count = polyCounts[f];
    if (!faceTags.empty()) {
      const std::string material = materialOf(f);
      if (material != lastMaterial) {
        ofs << "usemtl " << material << "\n";
        lastMaterial = material;
      }
    }
    ofs << "f";
    for (int i = 0; i < count; i++)
      ofs << " " << (polyConnects[offset + i] + 1);
    ofs << "\n";
    offset += count;
  }
  return bool(ofs);
}

ZSPACE_INLINE bool exportGraphToOBJ(const zObjGraph &graph,
                                    const std::string &filePath) {
  zFnGraph fnGraph(const_cast<zObjGraph &>(graph));
  zPointArray positions;
  zIntArray edgeConnects;
  fnGraph.getVertexPositions(positions);
  fnGraph.getEdgeData(edgeConnects);

  std::ofstream ofs(filePath);
  if (!ofs.is_open())
    return false;

  // Chain edges into polylines: open chains run between vertices whose
  // degree is not 2; whatever is left afterwards are closed loops.
  const int numVerts = int(positions.size());
  const int numEdges = int(edgeConnects.size() / 2);
  std::vector<std::vector<std::pair<int, int>>> adjacency(numVerts);
  for (int e = 0; e < numEdges; e++) {
    const int a = edgeConnects[2 * e], b = edgeConnects[2 * e + 1];
    adjacency[a].push_back({b, e});
    adjacency[b].push_back({a, e});
  }
  std::vector<bool> usedEdge(numEdges, false);
  std::vector<zIntArray> polylines;

  const auto walk = [&](int start, int firstEdge) {
    zIntArray polyline = {start};
    int current = start, edge = firstEdge;
    while (edge >= 0) {
      usedEdge[edge] = true;
      const int a = edgeConnects[2 * edge], b = edgeConnects[2 * edge + 1];
      current = (a == current) ? b : a;
      polyline.push_back(current);
      edge = -1;
      if (current == start || adjacency[current].size() != 2)
        break;
      for (const auto &nb : adjacency[current])
        if (!usedEdge[nb.second]) {
          edge = nb.second;
          break;
        }
    }
    polylines.push_back(polyline);
  };

  for (int v = 0; v < numVerts; v++)
    if (adjacency[v].size() != 2)
      for (const auto &nb : adjacency[v])
        if (!usedEdge[nb.second])
          walk(v, nb.second);
  for (int e = 0; e < numEdges; e++)
    if (!usedEdge[e])
      walk(edgeConnects[2 * e], e);

  ofs << "# zTsIO graph: " << numVerts << " vertices, " << numEdges
      << " edges as " << polylines.size() << " polylines\n";
  _writePositions(ofs, positions, zColorArray());
  for (const zIntArray &polyline : polylines) {
    ofs << "l";
    for (int v : polyline)
      ofs << " " << (v + 1);
    ofs << "\n";
  }
  return bool(ofs);
}

} // namespace zTsIO
} // namespace zSpace
