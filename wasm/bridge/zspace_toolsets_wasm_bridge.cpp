#include <emscripten/emscripten.h>
#include <zspace/toolsets.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdlib>
#include <limits>
#include <emscripten/heap.h>
#include "toolsets_binary_transfer.h"

namespace
{
using namespace zSpace;
using Json=nlohmann::json;
struct Session
{
    zTsCTubes tubes;
    std::unique_ptr<zTsVault> vault;
    std::unique_ptr<zTsPanelisation> panelisation;
    std::unique_ptr<zTsParametricStairs> stairs;
    Json stairsPostSmooth, stairsOptimisedGraph;
    std::string resultJson;
    Json resultData;
    std::vector<uint8_t> resultBuffer;
    std::vector<float> positions;
    std::vector<int> counts,connects;
};
std::map<int,std::unique_ptr<Session>> sessions;
int nextHandle=1;
std::string lastError;
Json* typedInput=nullptr;
bool typedOutput=false;
struct TransferScope {
    Json* previousInput;bool previousOutput;
    TransferScope(Json* input=nullptr,bool output=false):previousInput(typedInput),previousOutput(typedOutput){typedInput=input;typedOutput=output;}
    ~TransferScope(){typedInput=previousInput;typedOutput=previousOutput;}
};
Json readInput(const char* path,const char* message){
    if(typedInput)return std::move(*typedInput);
    std::ifstream stream(path?path:"");if(!stream)throw std::invalid_argument(message);
    Json input;stream>>input;return input;
}
const char* publish(Session& s,const Json& data){
    if(typedOutput){s.resultBuffer=zspaceTransfer::encode(data);return reinterpret_cast<const char*>(s.resultBuffer.data());}
    s.resultJson=data.dump();return s.resultJson.c_str();
}
template<class F>int invokeBuffer(const void* bytes,size_t size,F operation){
    int success=0;
    try {
        lastError.clear();auto address=reinterpret_cast<uintptr_t>(bytes);auto heap=emscripten_get_heap_size();
        if(!address||address>heap||size>heap-address)throw std::invalid_argument("Toolsets input buffer is outside memory.");
        auto input=zspaceTransfer::decode(bytes,size);TransferScope scope(&input);success=operation();
    }catch(const std::exception& e){lastError=e.what();}catch(...){lastError="Unknown toolsets transfer error.";}
    return success;
}
Session& session(int handle)
{
    auto it=sessions.find(handle);
    if(it==sessions.end()) throw std::invalid_argument("Invalid toolset handle.");
    return *it->second;
}
zTsVault& vault(int h) {
    auto& s=session(h);
    if(!s.vault) throw std::invalid_argument("Handle is not a Vault session.");
    return *s.vault;
}
template<class F> int guard(F function)
{
    try { lastError.clear(); function(); return 1; }
    catch(const std::exception& e) { lastError=e.what(); return 0; }
    catch(...) { lastError="Unknown toolset error."; return 0; }
}
zObjectGraph graph(const Json& input)
{
    auto positions=input.at("vertexPositions").get<std::vector<double>>();
    auto edges=input.at("edgeConnects").get<std::vector<int>>();
    if(positions.size()%3 || positions.size()>300000 || edges.size()%2) throw std::invalid_argument("Invalid graph arrays.");
    zPointArray p; for(size_t i=0;i<positions.size();i+=3)
    {
        for(int j=0;j<3;++j) if(!std::isfinite(positions[i+j]) || std::abs(positions[i+j])>1e20) throw std::invalid_argument("Invalid graph coordinate.");
        p.emplace_back(float(positions[i]),float(positions[i+1]),float(positions[i+2]));
    }
    for(int e:edges) if(e<0 || e>=int(p.size())) throw std::invalid_argument("Graph index out of range.");
    zObjectGraph result; zFnGraph(result).create(p,edges); return result;
}
}
extern "C"
{
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_abi_version() { return 1; }
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_transfer_version() { return 1; }
EMSCRIPTEN_KEEPALIVE void* zspace_alloc(size_t size) { return size?std::malloc(size):nullptr; }
EMSCRIPTEN_KEEPALIVE void zspace_free(void* pointer) { std::free(pointer); }
EMSCRIPTEN_KEEPALIVE int zspace_panelisation_create() {
    int handle=0;guard([&]{auto s=std::make_unique<Session>();s->panelisation=std::make_unique<zTsPanelisation>();handle=nextHandle++;sessions.emplace(handle,std::move(s));});return handle;
}
EMSCRIPTEN_KEEPALIVE int zspace_panelisation_compute(int handle,const char* path) {
    return guard([&]{
        auto& current=session(handle);if(!current.panelisation)throw std::invalid_argument("Handle is not a Panelisation session.");
        Json input=readInput(path,"Cannot open panelisation input.");
        std::vector<zPanelisationSource> sources;size_t total=0;
        for(const auto& item:input.at("meshes")) {
            auto values=item.at("vertexPositions").get<std::vector<double>>();auto counts=item.at("faceCounts").get<zIntArray>();auto connects=item.at("faceConnects").get<zIntArray>();
            total+=values.size()/3;
            if(values.empty()||values.size()%3||total>500000||connects.size()>2000000)throw std::invalid_argument("Invalid panel mesh size (limit 500,000 vertices).");
            size_t sum=0;for(int c:counts){if(c<3||c>4)throw std::invalid_argument("Panel inputs must be triangles or quads.");sum+=c;}
            if(sum!=connects.size())throw std::invalid_argument("Invalid panel connectivity length.");
            for(int v:connects)if(v<0||size_t(v)>=values.size()/3)throw std::invalid_argument("Invalid panel vertex index.");
            zPointArray p;for(size_t i=0;i<values.size();i+=3){for(int j=0;j<3;++j)if(!std::isfinite(values[i+j])||std::abs(values[i+j])>1e8)throw std::invalid_argument("Invalid panel coordinate.");p.emplace_back(float(values[i]),float(values[i+1]),float(values[i+2]));}
            sources.emplace_back();zFnMesh(sources.back().mesh).create(p,counts,connects);sources.back().creaseEdges=item.value("creaseEdges",zIntArray{});sources.back().panelSeamSeeds=item.value("panelSeamSeeds",zIntArray{});sources.back().topologyGuides=item.value("topologyGuides",false);sources.back().layoutMode=item.value("layoutMode",0);sources.back().alignmentSource=item.value("alignmentSource",-1);sources.back().alignmentDirection=item.value("alignmentDirection",std::vector<double>{});sources.back().staggerFraction=item.value("staggerFraction",0.5);sources.back().panelFrame=item.value("panelFrame",std::vector<double>{});
        }
        auto options=input.value("settings",Json::object());zPanelisationSettings settings;
        settings.width=options.value("width",1.0);settings.spacing=options.value("spacing",1.0);settings.offset=options.value("offset",0.0);settings.seedBand=options.value("seedBand",0.05);
        settings.weldTolerance=options.value("weldTolerance",0.0001);settings.creaseAngleDegrees=options.value("creaseAngle",35.0);settings.smallPanelArea=options.value("smallPanelArea",0.01);
        settings.axis=options.value("axis",0);settings.reverse=options.value("reverse",false);settings.useCreases=options.value("useCreases",true);
        auto candidate=std::make_unique<zTsPanelisation>();candidate->compute(sources,settings);
        Json panels=Json::array();for(const auto& p:candidate->panels())panels.push_back({{"source",p.source},{"station",p.station},{"area",p.area},{"column",p.column},{"chart",p.chart}});
        Json result={{"panels",panels},{"facePanels",candidate->facePanels()},{"faceSources",candidate->faceSources()},
            {"jointEdges",candidate->jointEdges()},{"jointKinds",candidate->jointKinds()},{"components",candidate->components()},
            {"weldedEdges",candidate->weldedEdges()},{"smallPanels",candidate->smallPanels()},{"axis",candidate->resolvedAxis()}};
        const auto& g=candidate->guides();result["guides"]={{"positions",g.positions},{"kinds",g.kinds},{"sources",g.sources},{"seeds",g.seeds},{"branches",g.branches},{"unsupportedCharts",g.unsupportedCharts},{"maxSampledWidth",g.maxSampledWidth}};
        current.panelisation=std::move(candidate);current.resultData=std::move(result);
    });
}
EMSCRIPTEN_KEEPALIVE const char* zspace_panelisation_result(int handle) {
    const char* value=nullptr;guard([&]{auto& s=session(handle);if(!s.panelisation)throw std::invalid_argument("Handle is not a Panelisation session.");value=publish(s,s.resultData);});return value;
}
EMSCRIPTEN_KEEPALIVE int zspace_panelisation_select_mesh(int handle) {
    return guard([&]{auto& s=session(handle);if(!s.panelisation)throw std::invalid_argument("Handle is not a Panelisation session.");
        zFnMesh fn(s.panelisation->mesh());zPointArray p;fn.getVertexPositions(p);s.positions.clear();for(auto v:p)s.positions.insert(s.positions.end(),{v.x,v.y,v.z});fn.getPolygonData(s.connects,s.counts);
    });
}
EMSCRIPTEN_KEEPALIVE const char* zspace_toolsets_last_error() { return lastError.c_str(); }
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_create()
{
    int handle=0;
    guard([&] { auto s=std::make_unique<Session>(); handle=nextHandle++; sessions.emplace(handle,std::move(s)); });
    return handle;
}
zParametricStairsSettings stairsSettings(const Json& input) {
    auto options=input.value("settings",Json::object());
    zParametricStairsSettings settings;settings.width=options.value("width",1.0);
        settings.numDivs=options.value("numDivs",1);settings.numSmoothDivs=options.value("numSmoothDivs",1);
        if(!options.contains("numDivs")&&options.contains("divisions"))settings.divisions=options.at("divisions").get<int>();
        settings.landingVertices=options.value("landingVertices",zIntArray{});settings.flightRiserCounts=options.value("flightRiserCounts",zIntArray{});settings.flightTreadDepths=options.value("flightTreadDepths",zDoubleArray{});
        settings.vertexWidths=options.value("vertexWidths",zDoubleArray{});settings.miterLimit=options.value("miterLimit",8.0);settings.vertexMass=options.value("vertexMass",1.0);
        settings.generateSteps=options.value("generateSteps",false);settings.treadDepth=options.value("treadDepth",.30);settings.riserHeight=options.value("riserHeight",.17);settings.landingDepth=options.value("landingDepth",1.0);
        settings.stepThickness=options.value("stepThickness",.08);settings.structureThickness=options.value("structureThickness",.15);settings.clearance=options.value("clearance",.05);
        if(input.contains("landingPolygons"))for(const auto& polygon:input.at("landingPolygons")){zPointArray p;for(size_t i=0;i<polygon.size();i+=3){if(i+2>=polygon.size())throw std::invalid_argument("Landing coordinates must be XYZ triples.");p.emplace_back(polygon[i].get<float>(),polygon[i+1].get<float>(),polygon[i+2].get<float>());}settings.landingPolygons.push_back(p);}
        settings.bestFit.qMin=options.value("qMin",0.01);settings.bestFit.qMax=options.value("qMax",1000.0);
        settings.bestFit.fixXY=options.value("fixXY",true);settings.bestFit.planConstraint=settings.bestFit.fixXY?zVaultPlanConstraint::FixXY:zVaultPlanConstraint::Free;
        settings.bestFit.geometryWeight=options.value("geometryWeight",0.5);settings.bestFit.maxIterations=options.value("maxIterations",2000);
    return settings;
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_create() {
    int h=0;guard([&]{auto s=std::make_unique<Session>();s->stairs=std::make_unique<zTsParametricStairs>();h=nextHandle++;sessions.emplace(h,std::move(s));});return h;
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_update(int h,const char* path) {
    return guard([&]{auto& current=session(h);if(!current.stairs)throw std::invalid_argument("Handle is not a Stairs session.");
        Json input=readInput(path,"Cannot open stairs input.");auto g=graph(input.at("graph"));auto options=input.value("settings",Json::object());
        auto settings=stairsSettings(input);
        auto candidate=std::make_unique<zTsParametricStairs>();candidate->update(g,settings);current.stairs=std::move(candidate);current.stairsPostSmooth=nullptr;
    });
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_optimise_graph(int h,const char* path) {
    return guard([&]{auto& current=session(h);if(!current.stairs)throw std::invalid_argument("Handle is not a Stairs session.");
        Json input=readInput(path,"Cannot open stairs graph optimisation input.");auto g=graph(input.at("graph"));auto settings=stairsSettings(input);auto range=input.value("ranges",Json::object());
        auto r=zTsParametricStairs::optimiseGraph(g,settings,range.value("treadMin",settings.treadDepth),range.value("treadMax",settings.treadDepth),range.value("riserMin",settings.riserHeight),range.value("riserMax",settings.riserHeight));
        zPointArray p;zIntArray e;zFnGraph fn(r.graph);fn.getVertexPositions(p);fn.getEdgeData(e);std::vector<float> xyz;for(auto v:p)xyz.insert(xyz.end(),{v.x,v.y,v.z});Json polygons=Json::array();for(auto& poly:r.landingPolygons){std::vector<float> out;for(auto v:poly)out.insert(out.end(),{v.x,v.y,v.z});polygons.push_back(out);}
        current.stairsOptimisedGraph={{"graph",{{"vertexPositions",xyz},{"edgeConnects",e}}},{"landingPolygons",polygons},{"landingOffsets",r.landingOffsets},{"bottomNodes",r.bottomNodes},{"landingVertices",r.landingVertices},{"adjustedLandings",r.adjustedLandings},{"flightRiserCounts",r.riserCounts},{"flightTreadDepths",r.treadDepths},{"flightActualRisers",r.riserHeights}};
    });
}
EMSCRIPTEN_KEEPALIVE const char* zspace_stairs_optimised_graph(int h){
    const char* result=nullptr;guard([&]{auto& s=session(h);if(!s.stairs||s.stairsOptimisedGraph.is_null())throw std::invalid_argument("Optimise a stairs graph first.");result=publish(s,s.stairsOptimisedGraph);});return result;
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_step(int h,int iterations) {
    return guard([&]{auto& s=session(h);if(!s.stairs)throw std::invalid_argument("Handle is not a Stairs session.");if(iterations<1||iterations>100)throw std::invalid_argument("Stairs step must be 1..100.");s.stairs->stepBestFit(iterations);s.stairsPostSmooth=nullptr;});
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_reset(int h) {
    return guard([&]{auto& s=session(h);if(!s.stairs)throw std::invalid_argument("Handle is not a Stairs session.");s.stairs->resetBestFit();s.stairsPostSmooth=nullptr;});
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_snap(int h,double geometryWeight) {
    return guard([&]{auto& s=session(h);if(!s.stairs)throw std::invalid_argument("Handle is not a Stairs session.");s.stairs->snapFixedDensities(geometryWeight);s.stairsPostSmooth=nullptr;});
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_previous(int h) {
    return guard([&]{auto& s=session(h);if(!s.stairs)throw std::invalid_argument("Handle is not a Stairs session.");s.stairs->previousBestFit();s.stairsPostSmooth=nullptr;});
}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_smooth_after(int h,int divisions) {
    return guard([&]{auto& s=session(h);if(!s.stairs)throw std::invalid_argument("Handle is not a Stairs session.");
        if(divisions<1||divisions>5)throw std::invalid_argument("Post-FDM smooth divisions must be 1..5.");
        s.stairs->smoothResult(divisions);s.stairsPostSmooth=nullptr;
    });
}
EMSCRIPTEN_KEEPALIVE const char* zspace_stairs_result(int h) {
    const char* result=nullptr;guard([&]{auto& s=session(h);if(!s.stairs)throw std::invalid_argument("Handle is not a Stairs session.");auto& t=*s.stairs;
        auto mesh=[](zObjectMesh& m){zFnMesh fn(m);zPointArray p;zIntArray c,f;fn.getVertexPositions(p);fn.getPolygonData(f,c);std::vector<float> xyz;for(auto v:p)xyz.insert(xyz.end(),{v.x,v.y,v.z});return Json{{"vertexPositions",xyz},{"faceCounts",c},{"faceConnects",f}};};
        auto r=t.report();Json data={{"lowPolyMesh",mesh(t.lowPolyMesh())},{"subdividedMesh",mesh(t.subdividedMesh())},
            {"centerMesh",mesh(t.centerMesh())},{"mesh",mesh(t.resultMesh())},{"supports",t.supports()},
            {"meshVertexGraphVertices",t.meshVertexGraphVertices()},{"centerGraphMeshVertices",t.centerGraphMeshVertices()},
            {"forceDensities",t.forceDensities()},{"status",r.status},{"iteration",r.iteration},{"converged",r.converged},
            {"residualNorm",r.residualNorm},{"maxResidual",r.maxResidual},{"rank",r.rank},{"activeEdges",r.activeEdges},{"freeVertices",r.freeVertices},
            {"rmsXYDisplacement",r.rmsXYDisplacement},{"maxXYDisplacement",r.maxXYDisplacement}};
        if(t.hasStairLayout()){
            data["stairMesh"]=mesh(t.stairMesh());data["stairFaceTypes"]=t.stairFaceTypes();data["flightRiserCounts"]=t.flightRiserCounts();data["flightActualRisers"]=t.flightActualRisers();data["minimumClearance"]=t.minimumClearance();
            zPointArray p;zIntArray e;zFnGraph fn(t.walkingGraph());fn.getVertexPositions(p);fn.getEdgeData(e);std::vector<float> xyz;for(auto v:p)xyz.insert(xyz.end(),{v.x,v.y,v.z});data["walkingGraph"]={{"vertexPositions",xyz},{"edgeConnects",e}};
        }
        if(!s.stairsPostSmooth.is_null())data["postSmoothMesh"]=s.stairsPostSmooth;
        data["edgeConnects"]=t.edgeConnects();data["activeEdgeIds"]=t.activeEdgeIds();
        auto vectors=[](const zVectorArray& values){std::vector<float> out;for(auto p:values)out.insert(out.end(),{p.x,p.y,p.z});return out;};
        data["axialForces"]=t.axialForces();data["residualForces"]=vectors(t.residualForces());data["supportReactions"]=vectors(t.supportReactions());
        result=publish(s,data);
    });return result;
}
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_destroy(int handle)
{ return guard([&] { session(handle); sessions.erase(handle); }); }
// File input avoids ccall's limited string stack for large graphs.
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_compute(int handle,const char* path)
{
    return guard([&] {
        auto& current=session(handle);
        if(current.vault || current.panelisation || current.stairs) throw std::invalid_argument("Handle is not a C-Tubes session.");
        Json input=readInput(path,"Cannot open C-Tubes input file."); auto candidate=std::make_unique<Session>();
        auto directrix=graph(input.at("directrix")); candidate->tubes.setDirectrix(directrix);
        if(input.contains("generatrices")) {
            zObjectGraphArray profiles;for(const auto& value:input.at("generatrices")) profiles.push_back(graph(value));
            candidate->tubes.setGeneratrices(profiles);
        } else { auto generatrix=graph(input.at("generatrix"));candidate->tubes.setGeneratrix(generatrix); }
        const Json profiles=input.contains("generatrices")?input.at("generatrices"):Json::array({input.at("generatrix")});
        for(int i=0;i<int(profiles.size());++i) if(profiles[i].contains("originVertex")) {
            int id=profiles[i].at("originVertex").get<int>(); const auto& positions=profiles[i].at("vertexPositions");
            if(id<0 || size_t(3*id+2)>=positions.size()) throw std::invalid_argument("Origin vertex is outside its generatrix.");
            candidate->tubes.setGeneratrixOrigin(i,zPoint(positions[3*id].get<float>(),positions[3*id+1].get<float>(),positions[3*id+2].get<float>()));
        }
        zCTubeSettings settings;
        settings.samples=input.value("samples",64); settings.spline=input.value("spline",true);
        settings.scale=input.value("scale",1.0); settings.twistRadians=input.value("twistRadians",0.0);
        settings.fitIterations=input.value("fitIterations",200); settings.fitTolerance=input.value("fitTolerance",0.001);
        settings.fitPlaneNormals=input.value("fitPlaneNormals",true);
        settings.targetSectionWeight=input.value("targetSectionWeight",1.0);
        settings.targetPlaneWeight=input.value("targetPlaneWeight",1.0);
        settings.apexSmoothnessWeight=input.value("apexSmoothnessWeight",0.0);
        settings.planeSmoothnessWeight=input.value("planeSmoothnessWeight",0.0);
        if(input.contains("sectionWeights")) settings.sectionWeights=input.at("sectionWeights").get<zDoubleArray>();
        if(input.contains("planeWeights")) settings.planeWeights=input.at("planeWeights").get<zDoubleArray>();
        settings.preserveOrientation=input.value("preserveOrientation",false);
        settings.apexScale=input.value("apexScale",1.0); candidate->tubes.setSettings(settings);
        int mode=input.value("mode",0); if(mode<0 || mode>1) throw std::invalid_argument("Unknown C-Tube mode.");
        if(input.contains("frames")) {
            if(mode!=0)throw std::invalid_argument("Fixed section frames use Generic C-Tube construction.");
            std::vector<zCTubeFrame> frames;
            for(const auto& f:input.at("frames")) {
                auto xyz=[&](const char* key) {const auto a=f.at(key).get<std::vector<double>>();
                    if(a.size()!=3)throw std::invalid_argument("Frame vectors require three coordinates.");
                    for(double v:a)if(!std::isfinite(v) || std::abs(v)>1e20)throw std::invalid_argument("Invalid frame coordinate.");
                    return zPoint(float(a[0]),float(a[1]),float(a[2]));};
                frames.push_back({xyz("origin"),xyz("x"),xyz("y"),xyz("normal")});
            }
            if(input.contains("sections")) {
                auto values=input.at("sections").get<std::vector<double>>();if(values.size()%3)throw std::invalid_argument("Invalid prescribed section coordinates.");
                zPointArray points;for(size_t i=0;i<values.size();i+=3){for(int k=0;k<3;++k)if(!std::isfinite(values[i+k]))throw std::invalid_argument("Invalid prescribed section coordinate.");points.emplace_back(float(values[i]),float(values[i+1]),float(values[i+2]));}
                candidate->tubes.computeTubeWithSections(frames,points);
            } else {
            zPointArray targets;
            if(input.contains("flangeTargets"))for(const auto& a:input.at("flangeTargets")) {
                auto p=a.get<std::vector<double>>();if(p.size()!=3)throw std::invalid_argument("Flange targets require three coordinates.");
                for(double v:p)if(!std::isfinite(v) || std::abs(v)>1e20)throw std::invalid_argument("Invalid flange target.");
                targets.emplace_back(float(p[0]),float(p[1]),float(p[2]));
            }
            candidate->tubes.computeTubeWithFrames(frames,targets,input.value("flangeVertex",-1));
            }
        } else candidate->tubes.computeTube(static_cast<zCTubeMode>(mode));
        sessions[handle]=std::move(candidate);
    });
}
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_unroll(int handle,double gap)
{ return guard([&] { auto& s=session(handle); if(s.vault || s.panelisation || s.stairs) throw std::invalid_argument("Handle is not a C-Tubes session."); auto& t=s.tubes; t.computeStrips(); t.unrollStrips(gap); }); }
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_branch_optimise(int handle,const char* path) {
    return guard([&]{auto& s=session(handle);if(s.vault||s.panelisation||s.stairs)throw std::invalid_argument("Not a C-Tubes session.");Json input=readInput(path,"Cannot open branching input.");
        auto xyz=[](const Json& a){auto p=a.get<std::vector<double>>();if(p.size()!=3)throw std::invalid_argument("Branch vectors need three coordinates.");for(double v:p)if(!std::isfinite(v)||std::abs(v)>1e8)throw std::invalid_argument("Invalid branch coordinate.");return zPoint(float(p[0]),float(p[1]),float(p[2]));};
        std::vector<zCTubeBranchInput> branches;for(const auto& b:input.at("branches")){branches.emplace_back();for(const auto& f:b.at("frames"))branches.back().frames.push_back({xyz(f.at("origin")),xyz(f.at("x")),xyz(f.at("y")),xyz(f.at("normal"))});for(const auto& p:b.at("centralTargets"))branches.back().centralTargets.push_back(xyz(p));}
        zCTubeBranchSettings settings;settings.outwardLength=input.value("outwardLength",.2);settings.stemLength=input.value("stemLength",.4);settings.iterations=input.value("iterations",400);settings.tolerance=input.value("tolerance",.001);settings.fixedDirectrices=input.value("fixedDirectrices",false);settings.fixXY=input.value("fixXY",false);
        zTsCTubeBranching optimiser;optimiser.optimise(branches,settings);Json output=Json::array();auto values=[](const zPoint& p){return Json::array({p.x,p.y,p.z});};
        for(const auto& b:optimiser.branches()){Json frames=Json::array(),sections=Json::array();for(const auto& f:b.frames)frames.push_back({{"origin",values(f.origin)},{"x",values(f.x)},{"y",values(f.y)},{"normal",values(f.normal)}});for(auto p:b.sections)for(float v:{p.x,p.y,p.z})sections.push_back(v);output.push_back({{"frames",frames},{"sections",sections}});}
        const auto& r=optimiser.report();s.resultData=Json({{"branches",output},{"report",{{"solver","L-BFGS-B"},{"gradient","Analytic chain rule"},{"iterations",r.iterations},{"rmsDisplacement",r.rmsDisplacement},{"maxDisplacement",r.maxDisplacement},{"objective",r.maxResidual},{"initialObjective",r.initialResidual},{"gradientCheckError",r.gradientCheckError},{"seamError",r.seamError},{"converged",r.converged},{"fixedDirectrices",settings.fixedDirectrices},{"fixXY",settings.fixXY}}}});
    });
}
EMSCRIPTEN_KEEPALIVE const char* zspace_ctubes_branch_result(int handle) {const char* result=nullptr;guard([&]{result=publish(session(handle),session(handle).resultData);});return result;}
EMSCRIPTEN_KEEPALIVE const char* zspace_ctubes_diagnostics(int handle)
{
    const char* result=nullptr;
    guard([&] {
        auto& s=session(handle);
        if(s.vault || s.panelisation || s.stairs) throw std::invalid_argument("Handle is not a C-Tubes session.");
        auto xyz=[](const zVector& p) { return Json::array({p.x,p.y,p.z}); };
        Json frames=Json::array(), areas=Json::array();
        for(const auto& f:s.tubes.frames()) frames.push_back({{"origin",xyz(f.origin)},{"x",xyz(f.x)},{"y",xyz(f.y)},{"normal",xyz(f.normal)}});
        for(const auto& a:s.tubes.stripAreas()) areas.push_back({{"sourceArea",a.sourceArea},{"unrolledArea",a.unrolledArea},{"difference",a.difference},{"deviationPercent",a.deviationPercent}});
        const auto& fit=s.tubes.fitReport();
        Json fitting={{"active",fit.active},{"converged",fit.converged},{"iterations",fit.iterations},{"initialRms",fit.initialRms},{"rms",fit.rms},{"maxError",fit.maxError},{"stations",fit.stations},{"apexScales",fit.apexScales}};
        fitting["initialObjective"]=fit.initialObjective; fitting["objective"]=fit.objective;
        fitting["solver"]="L-BFGS-B";fitting["gradientMethod"]="Analytic chain rule";fitting["gradientNorm"]=fit.gradientNorm;fitting["termination"]=fit.termination;fitting["evaluations"]=fit.evaluations;
        fitting["sectionEnergy"]=fit.sectionEnergy; fitting["planeEnergy"]=fit.planeEnergy;
        fitting["apexSmoothnessEnergy"]=fit.apexSmoothnessEnergy; fitting["planeSmoothnessEnergy"]=fit.planeSmoothnessEnergy;
        fitting["sectionRms"]=fit.sectionRms; fitting["planeAngleDegrees"]=fit.planeAngleDegrees;
        result=publish(s,Json({{"frames",frames},{"stripAreas",areas},{"fit",fitting}}));
    }); return result;
}
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_strip_count(int handle)
{
    int count=-1; guard([&] { auto& s=session(handle); if(s.vault || s.panelisation || s.stairs) throw std::invalid_argument("Handle is not a C-Tubes session."); count=int(s.tubes.stripMeshes().size()); }); return count;
}
// kind 0=tube, 1=3D strip, 2=unrolled strip. Buffers remain valid until the next mutation.
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_select_mesh(int handle,int kind,int index)
{
    return guard([&] {
        auto& s=session(handle); zObjectMesh* mesh=nullptr;
        if(s.vault || s.panelisation || s.stairs) throw std::invalid_argument("Handle is not a C-Tubes session.");
        if(kind==0) mesh=&s.tubes.tubeMesh();
        else if(kind==1) mesh=&s.tubes.stripMeshes().at(index);
        else if(kind==2) mesh=&s.tubes.unrolledMeshes().at(index);
        else throw std::invalid_argument("Unknown mesh kind.");
        zFnMesh fn(*mesh); zPointArray p; fn.getVertexPositions(p);
        s.positions.clear(); for(auto& v:p) s.positions.insert(s.positions.end(),{v.x,v.y,v.z});
        fn.getPolygonData(s.connects,s.counts);
    });
}
EMSCRIPTEN_KEEPALIVE const float* zspace_toolsets_positions(int h) { auto it=sessions.find(h); return it==sessions.end()?nullptr:it->second->positions.data(); }
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_positions_count(int h) { auto it=sessions.find(h); return it==sessions.end()?0:int(it->second->positions.size()); }
EMSCRIPTEN_KEEPALIVE const int* zspace_toolsets_face_counts(int h) { auto it=sessions.find(h); return it==sessions.end()?nullptr:it->second->counts.data(); }
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_face_count(int h) { auto it=sessions.find(h); return it==sessions.end()?0:int(it->second->counts.size()); }
EMSCRIPTEN_KEEPALIVE const int* zspace_toolsets_face_connects(int h) { auto it=sessions.find(h); return it==sessions.end()?nullptr:it->second->connects.data(); }
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_face_connects_count(int h) { auto it=sessions.find(h); return it==sessions.end()?0:int(it->second->connects.size()); }

EMSCRIPTEN_KEEPALIVE int zspace_vault_create() {
    int handle=0;guard([&]{auto s=std::make_unique<Session>();s->vault=std::make_unique<zTsVault>();handle=nextHandle++;sessions.emplace(handle,std::move(s));});return handle;
}
EMSCRIPTEN_KEEPALIVE int zspace_vault_fit(int h,const char* path) {
    return guard([&]{
        vault(h);Json input=readInput(path,"Cannot open Vault input.");auto candidate=std::make_unique<zTsVault>();
        auto data=input.at("mesh");auto coordinates=data.at("vertexPositions").get<std::vector<double>>();
        auto counts=data.at("faceCounts").get<zIntArray>();auto connects=data.at("faceConnects").get<zIntArray>();
        if(coordinates.empty()||coordinates.size()%3||coordinates.size()>30000||counts.size()>30000||connects.size()>120000) throw std::invalid_argument("Invalid mesh array sizes.");
        size_t total=0;for(int c:counts){if(c<3||c>10000)throw std::invalid_argument("Invalid polygon size.");total+=c;}
        if(total!=connects.size())throw std::invalid_argument("Polygon counts do not match connectivity.");
        zPointArray p;for(size_t i=0;i<coordinates.size();i+=3){
            for(int c=0;c<3;++c)if(!std::isfinite(coordinates[i+c])||std::abs(coordinates[i+c])>1e12)throw std::invalid_argument("Invalid coordinate.");
            p.emplace_back(float(coordinates[i]),float(coordinates[i+1]),float(coordinates[i+2]));
        }
        for(int v:connects)if(v<0||v>=int(p.size()))throw std::invalid_argument("Invalid vertex index.");
        zObjectMesh mesh;zFnMesh(mesh).create(p,counts,connects);candidate->setInputMesh(mesh);
        if(input.contains("referencePositions")){
            auto coords=input["referencePositions"].get<zDoubleArray>();
            if(coords.size()!=coordinates.size())throw std::invalid_argument("Reference positions must match input vertices.");
            zPointArray reference;for(size_t i=0;i<coords.size();i+=3)reference.emplace_back(float(coords[i]),float(coords[i+1]),float(coords[i+2]));
            candidate->setReferencePositions(reference);
        }
        candidate->setSupports(input.value("supports",zIntArray{}));
        // Resolve endpoint pairs against the toolsets edge list, never viewer edge IDs.
        auto edgeData=candidate->edgeConnects();zIntArray tension;
        for(auto pair:input.value("tensionEdges",Json::array())){
            if(!pair.is_array()||pair.size()!=2)throw std::invalid_argument("Expected tension endpoint pair.");
            int a=pair[0].get<int>(),b=pair[1].get<int>(),found=-1;
            for(int e=0;e<int(edgeData.size()/2);++e)if((edgeData[2*e]==a&&edgeData[2*e+1]==b)||(edgeData[2*e]==b&&edgeData[2*e+1]==a)){found=e;break;}
            if(found<0)throw std::invalid_argument("Tension edge does not exist in assigned mesh.");tension.push_back(found);
        }
        candidate->setTensionEdges(tension);
        if(input.contains("masses"))candidate->setVertexMasses(input["masses"].get<zDoubleArray>());
        else candidate->setVertexMass(input.value("mass",1.0));
        auto vector=[](const Json& v){if(!v.is_array()||v.size()!=3)throw std::invalid_argument("Expected xyz vector.");return zVector(v[0].get<float>(),v[1].get<float>(),v[2].get<float>());};
        candidate->setGravity(vector(input.value("gravity",Json::array({0,0,-1}))));
        if(input.contains("loads")){zVectorArray loads;for(auto& v:input["loads"])loads.push_back(vector(v));candidate->setVertexLoads(loads);}
        zVaultBestFitSettings settings;auto s=input.value("settings",Json::object());
        settings.qMin=s.value("qMin",0.01);settings.qMax=s.value("qMax",1000.0);settings.compressionOnly=s.value("compressionOnly",true);
        settings.fixXY=s.value("fixXY",false);
        if(s.contains("planConstraint")&&!s["planConstraint"].is_number_integer())throw std::invalid_argument("Plan constraint must be an integer mode.");
        settings.planConstraint=static_cast<zVaultPlanConstraint>(s.value("planConstraint",0));settings.xyPenaltyWeight=s.value("xyPenaltyWeight",1.0);
        settings.minMaxEdgeRatio=s.value("minMaxEdgeRatio",0.1);settings.edgeRegularizationStrength=s.value("edgeRegularizationStrength",0.0);
        settings.geometryWeight=s.value("geometryWeight",0.5);settings.absoluteTolerance=s.value("absoluteTolerance",1e-5);settings.relativeTolerance=s.value("relativeTolerance",0.01);
        settings.initialStep=s.value("initialStep",0.1);settings.integrationTolerance=s.value("integrationTolerance",1e-4);settings.maxIterations=s.value("maxIterations",2000);
        candidate->setBestFitSettings(settings);
        if(input.contains("qLower"))candidate->setForceDensityBounds(input["qLower"].get<zDoubleArray>(),input.at("qUpper").get<zDoubleArray>());
        candidate->computeBestFitForceDensities();session(h).vault=std::move(candidate);
    });
}
EMSCRIPTEN_KEEPALIVE int zspace_vault_step(int h,int iterations){return guard([&]{vault(h).stepBestFit(iterations);});}
EMSCRIPTEN_KEEPALIVE int zspace_vault_reset(int h){return guard([&]{vault(h).resetBestFit();});}
EMSCRIPTEN_KEEPALIVE int zspace_vault_previous(int h){return guard([&]{vault(h).previousBestFit();});}
EMSCRIPTEN_KEEPALIVE const char* zspace_vault_diagrams(int h,double tolerance){
    const char* result=nullptr;guard([&]{
        auto& v=vault(h);auto r=v.createForceFromForm({tolerance});
        auto meshJson=[](const zObjectMesh& input){
            zObjectMesh copy=input;zFnMesh fn(copy);zPointArray p;zIntArray c,n;
            fn.getVertexPositions(p);fn.getPolygonData(c,n);Json points=Json::array();
            for(auto x:p){points.push_back(x.x);points.push_back(x.y);points.push_back(x.z);}
            return Json{{"vertexPositions",points},{"faceCounts",n},{"faceConnects",c},{"edges",Json::array()},{"densities",Json::array()}};
        };
        Json form=meshJson(v.formMesh()),force=meshJson(v.forceMesh());
        zObjectMesh formObject=v.formMesh(),forceObject=v.forceMesh();
        zFnMesh source(v.resultMesh());
        for(zItMeshHalfEdge he(formObject);!he.end();he++){
            if(he.getId()>he.getSym().getId())continue;
            int a=he.getStartVertex().getId(),b=he.getVertex().getId(),edge=-1;
            if(!source.edgeExists(a,b,edge))throw std::runtime_error("Missing source edge.");
            double q=v.forceDensities()[edge];
            form["edges"].push_back(a);form["edges"].push_back(b);form["densities"].push_back(q);
            int mapped=v.formHalfEdgeToForceHalfEdge()[he.getId()];
            if(mapped>=0){zItMeshHalfEdge f(forceObject,mapped);
                force["edges"].push_back(f.getStartVertex().getId());force["edges"].push_back(f.getVertex().getId());force["densities"].push_back(q);
            }
        }
        Json data={{"form",form},{"force",force},{"status",r.status},{"exact",r.exact},
            {"horizontalResidual",r.horizontalResidual},{"closureError",r.closureError},{"maxEdgeMismatch",r.maxEdgeMismatch},
            {"formToForce",v.formHalfEdgeToForceHalfEdge()},{"forceToForm",v.forceHalfEdgeToFormHalfEdge()}};
        result=publish(session(h),std::move(data));
    });return result;
}
EMSCRIPTEN_KEEPALIVE int zspace_vault_select_mesh(int h){return guard([&]{
    auto& s=session(h);zFnMesh fn(vault(h).resultMesh());zPointArray p;fn.getVertexPositions(p);
    s.positions.clear();for(auto v:p)s.positions.insert(s.positions.end(),{v.x,v.y,v.z});fn.getPolygonData(s.connects,s.counts);
});}
EMSCRIPTEN_KEEPALIVE const char* zspace_vault_result(int h){
    const char* result=nullptr;guard([&]{auto& v=vault(h);auto& r=v.report();
        Json data={{"status",r.status},{"iteration",r.iteration},{"rank",r.rank},{"activeEdges",r.activeEdges},{"freeVertices",r.freeVertices},
            {"residualNorm",r.residualNorm},{"maxResidual",r.maxResidual},{"initialResidual",r.initialResidual},{"maxDisplacement",r.maxDisplacement},{"stepSize",r.stepSize},{"converged",r.converged},
            {"forceDensities",v.forceDensities()},{"axialForces",v.axialForces()},{"edgeConnects",v.edgeConnects()},{"activeEdgeIds",v.activeEdgeIds()},
            {"minMaxEdgeRatio",r.minMaxEdgeRatio},{"edgeRegularizationError",r.edgeRegularizationError},
            {"rmsXYDisplacement",r.rmsXYDisplacement},{"maxXYDisplacement",r.maxXYDisplacement}};
        auto vectors=[](const zVectorArray& values){Json out=Json::array();for(auto p:values){out.push_back(p.x);out.push_back(p.y);out.push_back(p.z);}return out;};
        data["residualForces"]=vectors(v.residualForces());data["supportReactions"]=vectors(v.supportReactions());
        result=publish(session(h),std::move(data));
    });return result;
}

// Additive binary APIs: legacy file/string entrypoints and ABI remain compatible.
EMSCRIPTEN_KEEPALIVE int zspace_toolsets_result_buffer_size(int h){auto it=sessions.find(h);return it==sessions.end()?0:int(it->second->resultBuffer.size());}
EMSCRIPTEN_KEEPALIVE int zspace_stairs_update_from_buffers(int h,const void* bytes,size_t size){return invokeBuffer(bytes,size,[&]{return zspace_stairs_update(h,nullptr);});}
EMSCRIPTEN_KEEPALIVE const void* zspace_stairs_result_buffers(int h){TransferScope scope(nullptr,true);return zspace_stairs_result(h);}
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_compute_from_buffers(int h,const void* bytes,size_t size){return invokeBuffer(bytes,size,[&]{return zspace_ctubes_compute(h,nullptr);});}
EMSCRIPTEN_KEEPALIVE int zspace_ctubes_branch_optimise_from_buffers(int h,const void* bytes,size_t size){return invokeBuffer(bytes,size,[&]{return zspace_ctubes_branch_optimise(h,nullptr);});}
EMSCRIPTEN_KEEPALIVE int zspace_vault_fit_from_buffers(int h,const void* bytes,size_t size){return invokeBuffer(bytes,size,[&]{return zspace_vault_fit(h,nullptr);});}
EMSCRIPTEN_KEEPALIVE int zspace_panelisation_compute_from_buffers(int h,const void* bytes,size_t size){return invokeBuffer(bytes,size,[&]{return zspace_panelisation_compute(h,nullptr);});}
EMSCRIPTEN_KEEPALIVE const void* zspace_ctubes_branch_result_buffers(int h){TransferScope scope(nullptr,true);return zspace_ctubes_branch_result(h);}
EMSCRIPTEN_KEEPALIVE const void* zspace_ctubes_diagnostics_buffers(int h){TransferScope scope(nullptr,true);return zspace_ctubes_diagnostics(h);}
EMSCRIPTEN_KEEPALIVE const void* zspace_vault_result_buffers(int h){TransferScope scope(nullptr,true);return zspace_vault_result(h);}
EMSCRIPTEN_KEEPALIVE const void* zspace_panelisation_result_buffers(int h){TransferScope scope(nullptr,true);return zspace_panelisation_result(h);}
EMSCRIPTEN_KEEPALIVE const void* zspace_vault_diagrams_buffers(int h,double tolerance){TransferScope scope(nullptr,true);return zspace_vault_diagrams(h,tolerance);}
}

extern "C" {
EMSCRIPTEN_KEEPALIVE int zspace_stairs_optimise_graph_from_buffers(int h,const void* bytes,size_t size){return invokeBuffer(bytes,size,[&]{return zspace_stairs_optimise_graph(h,nullptr);});}
EMSCRIPTEN_KEEPALIVE const void* zspace_stairs_optimised_graph_buffers(int h){TransferScope scope(nullptr,true);return zspace_stairs_optimised_graph(h);}
}
