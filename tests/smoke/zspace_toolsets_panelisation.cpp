#include <zspace/toolsets.h>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <fstream>
#include <nlohmann/json.hpp>
using namespace zSpace;
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(int argc,char** argv){try{
    if(argc>1){
        nlohmann::json input;std::ifstream stream(argv[1]);stream>>input;std::vector<zPanelisationSource> sources;
        for(const auto& m:input.at("meshes")){auto a=m.at("vertexPositions").get<std::vector<double>>();zPointArray p;for(size_t i=0;i<a.size();i+=3)p.emplace_back(float(a[i]),float(a[i+1]),float(a[i+2]));auto c=m.at("faceCounts").get<zIntArray>();auto f=m.at("faceConnects").get<zIntArray>();sources.emplace_back();zFnMesh(sources.back().mesh).create(p,c,f);sources.back().creaseEdges=m.at("creaseEdges").get<zIntArray>();}
        zPanelisationSettings settings;settings.weldTolerance=input.value("settings",nlohmann::json::object()).value("weldTolerance",0.0001);
        zTsPanelisation tool;tool.compute(sources,settings);double area=0;for(auto p:tool.panels())area+=p.area;
        std::cout<<"Real input: "<<tool.panels().size()<<" panels, "<<tool.components()<<" components, "<<tool.weldedEdges()<<" welded edges, "<<tool.smallPanels()<<" small panels, "<<area<<" m2\n";
        if(argc>2){zPointArray points;zIntArray c,f;zFnMesh(tool.mesh()).getVertexPositions(points);zFnMesh(tool.mesh()).getPolygonData(f,c);std::vector<double> xyz;for(auto v:points)xyz.insert(xyz.end(),{v.x,v.y,v.z});nlohmann::json panels=nlohmann::json::array();for(auto v:tool.panels())panels.push_back({{"source",v.source},{"chart",v.chart},{"station",v.station},{"column",v.column},{"area",v.area}});std::ofstream out(argv[2]);out<<nlohmann::json{{"vertexPositions",xyz},{"faceCounts",c},{"faceConnects",f},{"facePanels",tool.facePanels()},{"jointEdges",tool.jointEdges()},{"jointKinds",tool.jointKinds()},{"panels",panels}};}
        return 0;
    }
    std::vector<zPanelisationSource> sources(1);
    zPointArray p={{0,0,0},{0,1,0},{0,2,0},{2,0,0},{2,1,0},{2,2,0},{2,0,2},{2,1,2},{2,2,2}};
    zIntArray counts={4,4,4,4},connects={0,3,4,1,1,4,5,2,3,6,7,4,4,7,8,5};
    zFnMesh(sources[0].mesh).create(p,counts,connects);sources[0].creaseEdges={1,4,4,7};
    zTsPanelisation tool;zPanelisationSettings settings;settings.axis=1;tool.compute(sources,settings);
    check(tool.panels().size()==8,"Bent ribbon should have eight one-metre panels across two crease strips");
    double area=0;for(auto panel:tool.panels()){area+=panel.area;check(std::abs(panel.area-1)<1e-5,"Panel area on bent ribbon");}
    check(std::abs(area-8)<1e-5,"Area conservation");check(tool.components()==1,"One connected component");
    zPointArray out;zIntArray c,f;zFnMesh(tool.mesh()).getVertexPositions(out);zFnMesh(tool.mesh()).getPolygonData(f,c);
    bool bendCut=false;for(auto v:out)if(std::abs(v.x-2)<1e-5&&std::abs(v.z-1)<1e-5)bendCut=true;
    check(bendCut,"Distance follows bent surface beyond world-axis endpoint");
    std::map<std::pair<int,int>,int> edges;int cursor=0;
    for(int n:c){for(int j=0;j<n;++j)++edges[std::minmax(f[cursor+j],f[cursor+(j+1)%n])];cursor+=n;}
    for(auto item:edges){check(item.second<=2,"Manifold output");if(item.second==1){auto a=out[item.first.first],b=out[item.first.second];check((std::abs(a.y)<1e-5&&std::abs(b.y)<1e-5)||(std::abs(a.y-2)<1e-5&&std::abs(b.y-2)<1e-5)||(std::abs(a.x)<1e-5&&std::abs(b.x)<1e-5)||(std::abs(a.z-2)<1e-5&&std::abs(b.z-2)<1e-5),"No internal crack or T-junction");}}
    settings.offset=0.25;tool.compute(sources,settings);area=0;for(auto panel:tool.panels())area+=panel.area;check(std::abs(area-8)<1e-5,"Offset preserves complete surface coverage");
    settings.offset=0;settings.width=0.5;tool.compute(sources,settings);check(tool.panels().size()==16,"Y subdivides across the two crease strips");
    for(auto panel:tool.panels())check(std::abs(panel.area-0.5)<1e-5,"Independent X/Y panel dimensions");
    // Rotation cannot change an intrinsic layout (a wall need not face a world axis).
    for(auto& v:p){double x=v.x,y=v.y,z=v.z;v.x=float(0.36*x-0.8*y+0.48*z);v.y=float(0.48*x+0.6*y+0.64*z);v.z=float(-0.8*x+0.6*z);}
    zFnMesh(sources[0].mesh).create(p,counts,connects);tool.compute(sources,settings);check(tool.panels().size()==16,"Rotated surface retains X/Y grid");
    for(auto panel:tool.panels())check(std::abs(panel.area-0.5)<1e-5,"Rotated grid preserves panel areas");
    settings.width=0;bool widthRejected=false;try{tool.compute(sources,settings);}catch(const std::invalid_argument&){widthRejected=true;}check(widthRejected,"Reject zero panel width");settings.width=1;
    settings.spacing=0;bool rejected=false;try{tool.compute(sources,settings);}catch(const std::invalid_argument&){rejected=true;}check(rejected,"Reject zero spacing");
    std::cout<<"Panelisation: surface stationing, crease strips, area and continuity passed\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
