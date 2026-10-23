#include <emscripten/emscripten.h>
#include <emscripten/heap.h>
#include <zspace/zToolsets/statics/zTsForceCTubes.h>
#include "toolsets_binary_transfer.h"
#include <map>

namespace {
using namespace zSpace;
using Json=nlohmann::json;
std::map<int,std::vector<uint8_t>> results;
int nextHandle=2000000;std::string lastError;
template<class F> int guard(F fn){try{lastError.clear();fn();return 1;}catch(const std::exception& e){lastError=e.what();return 0;}}
std::vector<uint8_t>& session(int h){auto it=results.find(h);if(it==results.end())throw std::invalid_argument("Invalid Force C-Tubes handle.");return it->second;}
zObjectGraph graph(const Json& data){
    auto positions=data.at("vertexPositions").get<std::vector<double>>();
    if(positions.size()<6||positions.size()%3||positions.size()>384)throw std::invalid_argument("Expected 2..128 graph vertices.");
    zPointArray p;for(size_t i=0;i<positions.size();i+=3){for(int k=0;k<3;++k)if(!std::isfinite(positions[i+k])||std::abs(positions[i+k])>1e12)throw std::invalid_argument("Invalid graph position.");p.emplace_back(float(positions[i]),float(positions[i+1]),float(positions[i+2]));}
    const auto& raw=data.at("edgeConnects");if(!raw.is_array()||raw.size()%2||raw.size()>512)throw std::invalid_argument("Invalid edge buffer.");
    zIntArray e;for(const auto& x:raw){double i=x.get<double>();if(!std::isfinite(i)||std::floor(i)!=i||i<0||i>=double(p.size()))throw std::invalid_argument("Invalid graph index.");e.push_back(int(i));}
    zObjectGraph g;zFnGraph(g).create(p,e);return g;
}
Json xyz(zPoint p){return {p.x,p.y,p.z};}
Json compute(const Json& input){
    auto guide=graph(input.at("directrix")),profile=graph(input.at("generatrix"));const auto& force=input.at("force");auto form=graph(force.at("graph"));
    zTsForceCTubes tool;tool.setInput(guide,profile);zForceCTubeSettings s;
    const auto settings=input.value("settings",Json::object());
    s.material=settings.value("thinWall",true)?zForceCTubeMaterial::ThinWall:zForceCTubeMaterial::SolidPolygon;
    s.thickness=settings.value("thickness",s.thickness);s.forceScale=settings.value("forceScale",s.forceScale);s.referenceStress=settings.value("referenceStress",s.referenceStress);
    s.areaTolerance=settings.value("areaTolerance",s.areaTolerance);s.correspondingAreaTolerance=settings.value("correspondingAreaTolerance",s.correspondingAreaTolerance);
    s.smoothness=settings.value("smoothness",s.smoothness);s.minScale=settings.value("minScale",s.minScale);s.maxScale=settings.value("maxScale",s.maxScale);s.maxIterations=settings.value("maxIterations",s.maxIterations);tool.setSettings(s);
    tool.setForceData(form,force.at("memberFaceAreas").get<zDoubleArray>(),force.at("axialForces").get<zDoubleArray>(),force.at("pairedAreaDifferences").get<zDoubleArray>(),force.value("equilibrium",false));
    const auto baseline=tool.compute(false);const double inputArea=tool.sections().front().area;
    auto r=input.value("optimise",false)?tool.compute(true):baseline;
    const double areaRatio=tool.sections().front().area/inputArea;
    Json out={{"inputSectionArea",inputArea},{"seedScale",s.material==zForceCTubeMaterial::ThinWall?areaRatio:std::sqrt(areaRatio)},{"status",r.status},{"termination",r.termination},{"iterations",r.iterations},{"initialAreaError",r.initialAreaError},{"maxAreaError",r.maxAreaError},{"maxEccentricity",r.maxEccentricity},{"maxQuadWarp",r.maxQuadWarp},{"areaResolved",r.areaResolved},{"forceAreasResolved",r.forceAreasResolved},{"sourceEquilibrium",r.sourceEquilibrium},{"sectionChecksPassed",r.sectionChecksPassed},{"apexScales",tool.apexScales()},{"sections",Json::array()},{"frames",Json::array()},{"materialModel",s.material==zForceCTubeMaterial::ThinWall?"thin_wall_strips_joint_overlaps_neglected":"solid_polygon"}};
    zPointArray p;zIntArray counts,connects;zFnMesh fn(tool.tubeMesh());fn.getVertexPositions(p);fn.getPolygonData(connects,counts);std::vector<double> positions;for(auto v:p){positions.push_back(v.x);positions.push_back(v.y);positions.push_back(v.z);}out["mesh"]={{"vertexPositions",positions},{"faceCounts",counts},{"faceConnects",connects}};
    for(const auto& f:tool.frames())out["frames"].push_back({{"origin",xyz(f.origin)},{"x",xyz(f.x)},{"y",xyz(f.y)},{"normal",xyz(f.normal)}});
    for(const auto& c:tool.sections())out["sections"].push_back({{"area",c.area},{"targetArea",c.targetArea},{"relativeAreaError",c.relativeAreaError},{"centroid",xyz(c.centroid)},{"thrustPoint",xyz(c.thrustPoint)},{"eccentricity",c.eccentricity},{"ex",c.ex},{"ey",c.ey},{"kernUtilisation",c.kernUtilisation},{"minStress",c.minStress},{"maxStress",c.maxStress},{"stressUtilisation",c.stressUtilisation},{"intersectionValid",c.intersectionValid},{"compressionOnly",c.compressionOnly},{"withinStress",c.withinStress},{"status",c.status}});
    return out;
}
}
extern "C" {
EMSCRIPTEN_KEEPALIVE int zspace_forcectubes_create(){int h=nextHandle++;results.emplace(h,std::vector<uint8_t>{});return h;}
EMSCRIPTEN_KEEPALIVE int zspace_forcectubes_destroy(int h){return guard([&]{session(h);results.erase(h);});}
EMSCRIPTEN_KEEPALIVE const char* zspace_forcectubes_last_error(){return lastError.c_str();}
EMSCRIPTEN_KEEPALIVE int zspace_forcectubes_compute_from_buffers(int h,const void* bytes,size_t size){return guard([&]{
    auto& buffer=session(h);const size_t address=reinterpret_cast<uintptr_t>(bytes),heap=emscripten_get_heap_size();if(!address||address>heap||size>heap-address)throw std::invalid_argument("Force C-Tubes input outside WASM memory.");
    auto out=compute(zspaceTransfer::decode(bytes,size));auto next=zspaceTransfer::encode(out);buffer=std::move(next);
});}
EMSCRIPTEN_KEEPALIVE const void* zspace_forcectubes_result_buffers(int h){const void* p=nullptr;guard([&]{auto& b=session(h);if(b.empty())throw std::invalid_argument("Generate the Force C-Tube first.");p=b.data();});return p;}
EMSCRIPTEN_KEEPALIVE size_t zspace_forcectubes_result_buffer_size(int h){size_t n=0;guard([&]{n=session(h).size();});return n;}
}
