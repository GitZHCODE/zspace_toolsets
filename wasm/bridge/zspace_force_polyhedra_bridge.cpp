#include <emscripten/emscripten.h>
#include <emscripten/heap.h>
#include <zspace/zToolsets/statics/zTsForcePolyhedra.h>
#include "toolsets_binary_transfer.h"
#include <fstream>
#include <map>

namespace {
using namespace zSpace;
using Json=nlohmann::json;
struct Session {
    std::unique_ptr<zTsForcePolyhedra> tool;
    std::string json;
    std::vector<uint8_t> buffer;
    bool cells=false;
};
std::map<int,Session> sessions;
int nextHandle=1000000;
std::string lastError;
Session& session(int h){auto it=sessions.find(h);if(it==sessions.end())throw std::invalid_argument("Invalid Force Polyhedra handle.");return it->second;}
template<class F>int guard(F f){try{lastError.clear();f();return 1;}catch(const std::exception& e){lastError=e.what();return 0;}catch(...){lastError="Unknown force polyhedra error.";return 0;}}
std::vector<double> xyz(const zVectorArray& p){std::vector<double> out;for(const auto& v:p){out.push_back(v.x);out.push_back(v.y);out.push_back(v.z);}return out;}
zIntArray indices(const Json& data){
    if(!data.is_array()||data.size()>4096)throw std::invalid_argument("Invalid component index array.");
    zIntArray out;for(const auto& x:data){if(!x.is_number())throw std::invalid_argument("Component indices must be integers.");double n=x.get<double>();
        if(!std::isfinite(n)||std::floor(n)!=n||n<0||n>2147483647.0)throw std::invalid_argument("Invalid component index.");out.push_back(int(n));}return out;
}
zVectorArray vectors(const Json& data,size_t count){
    auto a=data.get<std::vector<double>>();if(a.size()!=3*count)throw std::invalid_argument("Expected XYZ triples.");
    zVectorArray out;for(size_t i=0;i<a.size();i+=3){for(int k=0;k<3;++k)if(!std::isfinite(a[i+k])||std::abs(a[i+k])>1e20)throw std::invalid_argument("Invalid vector coordinate.");out.emplace_back(float(a[i]),float(a[i+1]),float(a[i+2]));}return out;
}
std::unique_ptr<zTsForcePolyhedra> compute(const Json& input,bool& cells){
    const auto& graph=input.at("graph");auto a=graph.at("vertexPositions").get<std::vector<double>>();
    if(a.empty()||a.size()%3||a.size()>3072)throw std::invalid_argument("Expected 1..1024 graph vertices.");
    zPointArray p=vectors(graph.at("vertexPositions"),a.size()/3);auto edges=indices(graph.at("edgeConnects"));
    if(edges.size()%2||edges.size()>4096)throw std::invalid_argument("Invalid graph edge array.");
    for(int i:edges)if(i<0||i>=int(p.size()))throw std::invalid_argument("Graph edge index out of range.");
    zObjectGraph g;zFnGraph(g).create(p,edges);auto tool=std::make_unique<zTsForcePolyhedra>();tool->setInputGraph(g);
    if(input.contains("supports"))tool->setSupports(indices(input.at("supports")));
    if(input.contains("tensionEdges"))tool->setTensionEdges(indices(input.at("tensionEdges")));
    if(input.contains("edgeForceMagnitudes"))tool->setEdgeForceMagnitudes(input.at("edgeForceMagnitudes").get<zDoubleArray>());
    if(input.contains("masses"))tool->setVertexMasses(input.at("masses").get<zDoubleArray>());
    if(input.contains("mass"))tool->setVertexMass(input.at("mass").get<double>());
    if(input.contains("gravity"))tool->setGravity(vectors(input.at("gravity"),1)[0]);
    if(input.contains("loads"))tool->setVertexLoads(vectors(input.at("loads"),p.size()));
    zForcePolyhedraSettings settings;
    if(input.contains("settings")){
        auto s=input.at("settings");
        settings.equilibriumTolerance=s.value("equilibriumTolerance",settings.equilibriumTolerance);
        settings.areaTolerance=s.value("areaTolerance",settings.areaTolerance);
        settings.correspondingAreaTolerance=s.value("correspondingAreaTolerance",settings.correspondingAreaTolerance);
        settings.normalTolerance=s.value("normalTolerance",settings.normalTolerance);
        settings.forceScale=s.value("forceScale",settings.forceScale);
        settings.maxIterations=s.value("maxIterations",settings.maxIterations);
        settings.formWeight=s.value("formWeight",settings.formWeight);
        settings.angleTolerance=s.value("angleTolerance",settings.angleTolerance);
        settings.relaxation=s.value("relaxation",settings.relaxation);
        settings.minMaxEdgeRatio=s.value("minMaxEdgeRatio",settings.minMaxEdgeRatio);
        settings.areaIterationsPerStep=s.value("areaIterationsPerStep",settings.areaIterationsPerStep);
        settings.timeStep=s.value("timeStep",settings.timeStep);
        settings.maxVertexStepRatio=s.value("maxVertexStepRatio",settings.maxVertexStepRatio);
        settings.normalStepDegrees=s.value("normalStepDegrees",settings.normalStepDegrees);
    }
    tool->setSettings(settings);
    if(input.contains("forceDensities"))throw std::invalid_argument("Force-density form finding is not part of this 3DGS toolset.");
    tool->createForcePolyhedra();cells=true;
    if(input.value("solveEquilibrium",false))tool->solveEquilibrium(input.value("equilibriumIterations",2000));return tool;
}
Json result(Session& s){
    if(!s.tool)throw std::invalid_argument("Compute force polyhedra first.");auto& t=*s.tool;auto r=t.report();
    Json out={{"status",r.status},{"equilibrium",r.equilibrium},{"cellsConverged",r.cellsConverged},
        {"areaPursuitResolved",r.areaPursuitResolved},{"maxCorrespondingAreaDifference",r.maxCorrespondingAreaDifference},{"worstAreaDifferenceEdge",r.worstAreaDifferenceEdge},
        {"iterations",r.iterations},{"maxResidual",r.maxResidual},{"residualNorm",r.residualNorm},
        {"maxAreaError",r.maxAreaError},{"maxNormalError",r.maxNormalError},{"maxPlanarityError",r.maxPlanarityError},
        {"axialForces",t.axialForces()},{"equilibriumIterations",r.equilibriumIterations},
        {"maxAngleDeviation",r.maxAngleDeviation},{"maxContactAreaMismatch",r.maxContactAreaMismatch},
        {"integration","zFnMeshDynamics + zFnGraphDynamics"},{"maxCellVertexStep",r.maxCellVertexStep},{"maxFormVertexStep",r.maxFormVertexStep},
        {"worstAreaCell",r.worstAreaCell},{"worstAngleEdge",r.worstAngleEdge},{"worstMotionCell",r.worstMotionCell},
        {"residualForces",xyz(t.residualForces())},{"supportReactions",xyz(t.supportReactions())}};
    zObjectGraph graph=t.resultGraph();zPointArray p;zFnGraph(graph).getVertexPositions(p);
    out["graph"]={{"vertexPositions",xyz(p)},{"edgeConnects",t.edgeConnects()}};out["cells"]=Json::array();
    if(s.cells)for(int i=0;i<int(p.size());++i){
        const auto& cr=t.cellReport(i);zObjectMesh mesh=t.forceCell(i);zFnMesh fn(mesh);zPointArray vertices;zIntArray counts,connects;
        fn.getVertexPositions(vertices);fn.getPolygonData(connects,counts);
        Json cell={{"vertexId",i},{"status",cr.status},{"converged",cr.converged},{"forceRank",cr.forceRank},
            {"iterations",cr.iterations},{"closureError",cr.closureError},{"maxAreaError",cr.maxAreaError},
            {"maxNormalError",cr.maxNormalError},{"maxPlanarityError",cr.maxPlanarityError},
            {"maxVertexStep",cr.maxVertexStep},{"dynamicsSteps",cr.dynamicsSteps},
            {"vertexPositions",xyz(vertices)},{"faceCounts",counts},{"faceConnects",connects},{"faces",Json::array()}};
        for(const auto& f:t.cellFaces(i)){
            Json face={{"type",int(f.type)},{"normal",{f.normal.x,f.normal.y,f.normal.z}},
                {"targetArea",f.targetArea},{"area",f.area},{"sources",Json::array()}};
            for(auto source:f.sources)face["sources"].push_back({{"type",int(source.type)},{"id",source.id},{"magnitude",source.magnitude}});
            cell["faces"].push_back(std::move(face));
        }
        out["cells"].push_back(std::move(cell));
    }
    return out;
}
}
extern "C" {
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_create(){int h=0;guard([&]{h=nextHandle++;sessions.emplace(h,Session{});});return h;}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_destroy(int h){return guard([&]{session(h);sessions.erase(h);});}
EMSCRIPTEN_KEEPALIVE const char* zspace_forcepolyhedra_last_error(){return lastError.c_str();}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_compute(int h,const char* path){return guard([&]{
    auto& s=session(h);std::ifstream stream(path?path:"");if(!stream)throw std::invalid_argument("Cannot read force polyhedra input.");
    Json input;stream>>input;bool cells=false;auto tool=compute(input,cells);s.tool=std::move(tool);s.cells=cells;s.buffer.clear();s.json.clear();
});}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_compute_from_buffers(int h,const void* bytes,size_t length){return guard([&]{
    auto& s=session(h);size_t address=reinterpret_cast<uintptr_t>(bytes),heap=emscripten_get_heap_size();
    if(!address||address>heap||length>heap-address)throw std::invalid_argument("Force polyhedra input outside WASM memory.");
    auto input=zspaceTransfer::decode(bytes,length);bool cells=false;auto tool=compute(input,cells);s.tool=std::move(tool);s.cells=cells;s.buffer.clear();s.json.clear();
});}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_step(int h,int iterations){return guard([&]{auto& s=session(h);if(!s.tool)throw std::invalid_argument("Compute force polyhedra first.");s.tool->stepAreaPursuit(iterations);s.buffer.clear();s.json.clear();});}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_equilibrium(int h,int iterations){return guard([&]{auto& s=session(h);if(!s.tool)throw std::invalid_argument("Compute force polyhedra first.");s.tool->stepEquilibrium(iterations);s.buffer.clear();s.json.clear();});}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_set_form_weight(int h,double weight){return guard([&]{auto& s=session(h);if(!s.tool)throw std::invalid_argument("Compute force polyhedra first.");s.tool->setFormWeight(weight);s.buffer.clear();s.json.clear();});}
EMSCRIPTEN_KEEPALIVE int zspace_forcepolyhedra_set_form_from_buffers(int h,const void* bytes,size_t length){return guard([&]{
    auto& s=session(h);if(!s.tool)throw std::invalid_argument("Compute force polyhedra first.");size_t address=reinterpret_cast<uintptr_t>(bytes),heap=emscripten_get_heap_size();
    if(!address||address>heap||length>heap-address)throw std::invalid_argument("Form input outside WASM memory.");
    auto input=zspaceTransfer::decode(bytes,length);zObjectGraph graph=s.tool->resultGraph();s.tool->setFormPositions(vectors(input.at("vertexPositions"),zFnGraph(graph).numVertices()));s.buffer.clear();s.json.clear();
});}
EMSCRIPTEN_KEEPALIVE const char* zspace_forcepolyhedra_result(int h){const char* p=nullptr;guard([&]{auto& s=session(h);s.json=result(s).dump();p=s.json.c_str();});return p;}
EMSCRIPTEN_KEEPALIVE const void* zspace_forcepolyhedra_result_buffers(int h){const void* p=nullptr;guard([&]{auto& s=session(h);s.buffer=zspaceTransfer::encode(result(s));p=s.buffer.data();});return p;}
EMSCRIPTEN_KEEPALIVE size_t zspace_forcepolyhedra_result_buffer_size(int h){size_t size=0;guard([&]{size=session(h).buffer.size();});return size;}
}
