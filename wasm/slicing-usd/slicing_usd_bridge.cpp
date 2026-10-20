#include <emscripten/emscripten.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>
#include "toolsets_binary_transfer.h"
#include <pxr/usd/usd/stage.h>
#include <pxr/usd/usd/primRange.h>
#include <pxr/usd/usdGeom/mesh.h>
#include <pxr/usd/usdGeom/xformCache.h>
#include <pxr/usd/usdGeom/primvarsAPI.h>
#include <pxr/usd/usdGeom/metrics.h>

namespace {
using namespace zSpace;
using Json = nlohmann::json;
std::string error;
std::vector<uint8_t> bytes;
Json flatten(const zPointArray& points) {
  Json out = Json::array();
  for (const auto& p : points) {out.push_back(p.x); out.push_back(p.y); out.push_back(p.z);}
  return out;
}
Json mesh(zObjMesh& object) {
  zPointArray positions; zIntArray counts, connects;
  zFnMesh fn(object); fn.getVertexPositions(positions); fn.getPolygonData(connects, counts);
  return {{"vertexPositions",flatten(positions)},{"faceCounts",counts},{"faceConnects",connects}};
}
Json read(const char* path) {
  auto stage = pxr::UsdStage::Open(path);
  if (!stage) throw std::runtime_error("Cannot open the USD stage. Use a self-contained file or USDZ package.");
  const auto block = stage->GetPrimAtPath(pxr::SdfPath("/block"));
  if (block && (pxr::UsdGeomPrimvarsAPI(block).GetPrimvar(pxr::TfToken("FaceTag")) ||
      block.GetAttribute(pxr::TfToken("FaceTag")) || stage->GetRootLayer()->GetCustomLayerData().count("slicer"))) {
    zTsIO::SlicerUSDInput input; std::string diagnostic;
    if (!zTsIO::readSlicerInputUSD(path,input,&diagnostic)) throw std::runtime_error(diagnostic);
    const auto original=Json::parse(input.SettingsJSON);
    const auto& s=original.at("slicer"); const auto& p=original.at("printSyn");
    Json settings=Json::object();
    const std::pair<const char*,const char*> slicerKeys[] = {
      {"LayerSpacing","layerSpacing"},{"BoundaryWidth","boundaryWidths"},{"FieldCellSize","fieldCellSize"},
      {"DoubleBoundary","doubleBoundary"},{"OverlapWidth","overlapWidth"},{"FieldScaleUp","fieldScaleUp"},
      {"BoundarySlotOffset","boundarySlotOffset"},{"BracingSlotEven","bracingSlotEven"},{"BracingSlotOdd","bracingSlotOdd"},
      {"HTraceStepFraction","traceStepFraction"},{"HTraceMaxSteps","traceMaxSteps"},{"HEqualiseIterations","equaliseIterations"},
      {"HArapIterations","arapIterations"},{"HUnrollStart","unrollStart"}
    };
    for(const auto& entry:slicerKeys) if(s.contains(entry.first))settings[entry.second]=s[entry.first];
    const std::pair<const char*,const char*> printKeys[] = {{"ResampleLength","resampleLength"},{"FeatureAngleTol","featureAngle"},{"SeamOffset","seamOffset"},{"SeamSize","seamSize"}};
    for(const auto& entry:printKeys) if(p.contains(entry.first))settings[entry.second]=p[entry.first];
    settings["bracingWidth"]=s.at("BracingWidth").at(0);
    settings["bracingCurveWidths"]=input.BracingPrintWidths;
    zPointArray positions; zIntArray edges;
    for(auto& curve:input.BracingCurves) {
      zPointArray cp; zIntArray ce; zFnGraph fn(curve); fn.getVertexPositions(cp); fn.getEdgeData(ce);
      const int offset=positions.size(); positions.insert(positions.end(),cp.begin(),cp.end());
      for(int id:ce)edges.push_back(id+offset);
    }
    return {{"mesh",mesh(input.Block)},{"faceTags",input.FaceTags},{"faceWidthIds",input.FaceWidthIds},
      {"anchor",{input.StartPoint.x,input.StartPoint.y,input.StartPoint.z}},
      {"bracing",{{"vertexPositions",flatten(positions)},{"edgeConnects",edges}}},{"settings",settings},
      {"metersPerUnit",input.MetersPerUnit},{"upAxis",input.UpAxis}};
  }
  std::vector<pxr::UsdPrim> meshes;
  for(const auto& prim:stage->Traverse()) if(prim.IsA<pxr::UsdGeomMesh>())meshes.push_back(prim);
  pxr::UsdPrim selected;
  if(block && block.IsA<pxr::UsdGeomMesh>())selected=block;
  else if(meshes.size()==1)selected=meshes.front();
  else throw std::runtime_error("USD must contain one mesh, or name the slicing mesh /block.");
  pxr::UsdGeomMesh object(selected); pxr::UsdGeomXformCache cache;
  pxr::VtArray<pxr::GfVec3f> points; pxr::VtArray<int> counts,connects,holes;
  if(!object.GetPointsAttr().Get(&points) || !object.GetFaceVertexCountsAttr().Get(&counts) || !object.GetFaceVertexIndicesAttr().Get(&connects))throw std::runtime_error("USD mesh has no polygon geometry.");
  object.GetHoleIndicesAttr().Get(&holes); if(!holes.empty())throw std::runtime_error("Bake USD hole faces before slicing.");
  const auto transform=cache.GetLocalToWorldTransform(selected);
  Json positions=Json::array();
  for(const auto& point:points) {const auto v=transform.Transform(pxr::GfVec3d(point));if(pxr::UsdGeomGetStageUpAxis(stage)==pxr::UsdGeomTokens->y){positions.push_back(v[0]);positions.push_back(-v[2]);positions.push_back(v[1]);}else{for(int k=0;k<3;++k)positions.push_back(v[k]);}}
  zIntArray indices(connects.begin(),connects.end());
  pxr::TfToken orientation;object.GetOrientationAttr().Get(&orientation);
  if((orientation==pxr::UsdGeomTokens->leftHanded)!=(transform.GetDeterminant()<0)) {
    size_t offset=0;for(int count:counts){if(count<3||offset+count>indices.size())throw std::runtime_error("Invalid USD polygon topology.");std::reverse(indices.begin()+offset,indices.begin()+offset+count);offset+=count;}
  }
  return {{"mesh",{{"vertexPositions",positions},{"faceCounts",zIntArray(counts.begin(),counts.end())},{"faceConnects",indices}}},
    {"metersPerUnit",pxr::UsdGeomGetStageMetersPerUnit(stage)},{"upAxis",pxr::UsdGeomGetStageUpAxis(stage).GetString()}};
}
}
extern "C" {
EMSCRIPTEN_KEEPALIVE int zspace_slicing_usd_read(const char* path) {
  try {error.clear(); auto next=zspaceTransfer::encode(read(path));bytes=std::move(next);return 1;}
  catch(const std::exception& e){error=e.what();return 0;}
}
EMSCRIPTEN_KEEPALIVE const uint8_t* zspace_slicing_usd_result(){return bytes.data();}
EMSCRIPTEN_KEEPALIVE size_t zspace_slicing_usd_result_size(){return bytes.size();}
EMSCRIPTEN_KEEPALIVE const char* zspace_slicing_usd_error(){return error.c_str();}
}
