#include <zspace/toolsets.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace zSpace;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
zObjectGraph graph(zPointArray p,zIntArray e){zObjectGraph g;zFnGraph(g).create(p,e);return g;}
int main(){try{
    auto guide=graph({{0,0,0},{0,0,2},{0,0,4}},{0,1,1,2});
    auto rectangle=graph({{-1,-2,0},{1,-2,0},{1,2,0},{-1,2,0}},{2,3,0,1,3,0,1,2});
    auto thrust=graph({{1.0f/3,0,0},{1.0f/3,0,2},{1.0f/3,0,4}},{0,1,1,2});
    zTsForceCTubes tool;tool.setInput(guide,rectangle);zForceCTubeSettings s;s.material=zForceCTubeMaterial::SolidPolygon;s.smoothness=0;tool.setSettings(s);
    tool.setForceData(thrust,{8,8},{-8,-8},{0,0},true);auto report=tool.compute();
    require(report.areaResolved&&report.forceAreasResolved,"Exact rectangle area fit");
    for(const auto& c:tool.sections()){
        require(c.intersectionValid&&std::abs(c.area-8)<1e-5,"Analytical rectangle area");
        require(std::abs(c.eccentricity-1.0/3)<1e-6,"Section-plane eccentricity");
        require(std::abs(c.kernUtilisation-1)<1e-5&&std::abs(c.minStress)<1e-5&&std::abs(c.maxStress-2)<1e-5,"Rectangle middle-third boundary and stress");
        require(!c.withinStress,"Eccentricity and stress checks are independent");
    }
    // Actual thin-wall material area is perimeter * thickness, not enclosed area.
    s.material=zForceCTubeMaterial::ThinWall;s.thickness=.1;tool.setSettings(s);tool.compute();
    require(std::abs(tool.sections()[0].area-1.2)<1e-5,"Thin-wall rectangle strip area");
    auto tee=graph({{0,0,0},{2,0,0},{0,2,0},{-1,0,0}},{0,1,0,2,0,3});
    zTsForceCTubes branched;
    branched.setInput(guide,tee);branched.setSettings(s);branched.setForceData(guide,{.5,.5},{-.5,-.5},{0,0},true);branched.compute();
    for(const auto& c:branched.sections())require(std::abs(c.area-.5)<1e-5,"Branched section uses material edges once despite split mesh seams");
    s.material=zForceCTubeMaterial::SolidPolygon;tool.setSettings(s);tool.setForceData(guide,{2,8},{-2,-8},{0,.2},false);
    auto generated=tool.compute();auto fitted=tool.compute(true);
    require(fitted.maxAreaError<.01&&fitted.maxAreaError<generated.maxAreaError*.01,"Area optimisation improves nonuniform force targets");
    require(!fitted.forceAreasResolved&&!fitted.sourceEquilibrium&&fitted.status=="provisional_force_snapshot","Unresolved force snapshot cannot claim solved");
    require(fitted.maxQuadWarp<1e-5,"Projection construction retains planar strip quads");
    require(tool.apexScales().size()==2,"Per-segment apex controls");
    // Oblique Rhino-style section anchored on a non-centroid vertex. Neither its
    // world plane, roll, nor anchor should be replaced with an automatic frame.
    auto anchored=graph({{0,0,0},{2,0,-1},{2,3,-1},{0,3,0}},{0,1,1,2,2,3,3,0});
    zTsForceCTubes oriented;oriented.setInput(guide,anchored);oriented.setSettings(s);oriented.setForceData(guide,{8,8},{-8,-8},{0,0},true);oriented.compute();
    zPointArray actual;zFnMesh(oriented.tubeMesh()).getVertexPositions(actual);
    zPointArray expected;zFnGraph(anchored).getVertexPositions(expected);
    for(int k=0;k<4;++k)require((actual[k]-expected[k]).length()<1e-5,"Preserve imported section orientation and anchor");
    auto missing=graph({{10,0,10},{10,0,12},{10,0,14}},{0,1,1,2});tool.setForceData(missing,{8,8},{-8,-8},{0,0},true);tool.compute();
    require(!tool.sections()[0].intersectionValid&&!tool.report().sectionChecksPassed,"Missing intersection stays unresolved");
    bool rejected=false;try{tool.setForceData(guide,{0,8},{-2,-8},{0,0},true);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"Reject zero force-face target");
    rejected=false;s.thickness=-1;try{tool.setSettings(s);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"Reject invalid material thickness");
    std::cout<<"Force C-Tubes analytical sections, projection sizing and failure states passed.\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
