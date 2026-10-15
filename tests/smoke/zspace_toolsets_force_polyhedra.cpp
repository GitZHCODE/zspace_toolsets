#include <zspace/toolsets.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
using namespace zSpace;
void check(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,"Expected invalid input rejection");}
zObjectGraph makeGraph(zPointArray p,zIntArray e){zObjectGraph g;zFnGraph(g).create(p,e);return g;}
zPointArray points(const zTsForcePolyhedra& t){auto g=t.resultGraph();zPointArray p;zFnGraph(g).getVertexPositions(p);return p;}
void print(const zForcePolyhedraReport& r){std::cout<<r.status<<" residual "<<r.maxResidual<<" angle "<<r.maxAngleDeviation<<" area "<<r.maxAreaError<<" normal "<<r.maxNormalError<<" plane "<<r.maxPlanarityError<<" contact "<<r.maxContactAreaMismatch<<" steps "<<r.equilibriumIterations<<'\n';}
void inspect(zTsForcePolyhedra& tool,int n){
    for(int i=0;i<n;++i){
        auto mesh=tool.forceCell(i);zFnMesh fn(mesh);zPointArray p;fn.getVertexPositions(p);zIntArray counts,connects;fn.getPolygonData(connects,counts);
        for(auto x:p)check(std::isfinite(x.x)&&std::isfinite(x.y)&&std::isfinite(x.z),"Finite output cell");
        check(counts.size()==tool.cellFaces(i).size(),"Export retains face-to-source numbering");
        std::map<std::pair<int,int>,int> incidence,orientation;
        size_t offset=0;for(size_t j=0;j<counts.size();++j){
            zVector av(0,0,0);for(int k=0;k<counts[j];++k){int u=connects[offset+k],v=connects[offset+(k+1)%counts[j]];auto a=p[u],b=p[v];av+=a^b;auto edge=std::make_pair(std::min(u,v),std::max(u,v));++incidence[edge];orientation[edge]+=u<v?1:-1;}av*=.5f;
            check(std::abs(av.length()-tool.cellFaces(i)[j].area)<1e-4*std::max(1.0,tool.cellFaces(i)[j].area),"Exported polygon area matches solver");offset+=counts[j];
        }
        for(auto edge:incidence)check(edge.second==2&&orientation[edge.first]==0,"Force-cell mesh is closed with consistent winding");
    }
}
int main(){try{
    auto arch=makeGraph({{-1,0,0},{0,0,1},{1,0,0}},{0,1,1,2});
    zTsForcePolyhedra tool;tool.setInputGraph(arch);tool.setSupports({0,2});tool.setVertexMasses({1,2,1});tool.setEdgeForceMagnitudes({std::sqrt(2.0),std::sqrt(2.0)});
    auto r=tool.createForcePolyhedra();print(r);check(r.equilibrium,"Analytic arch equilibrium and reciprocity");inspect(tool,3);
    check(std::abs(tool.supportReactions()[0].z-2)<1e-6,"Supports carry their own self weight");
    check(tool.cellReport(1).forceRank==2,"Planar node rank");int virtuals=0;for(auto f:tool.cellFaces(1))if(f.type==zForcePolyhedraFaceType::Virtual)++virtuals;check(virtuals==2,"Two virtual caps");
    auto cable=makeGraph({{-1,0,0},{0,0,-1},{1,0,0}},{0,1,1,2});tool.setInputGraph(cable);tool.setSupports({0,2});tool.setVertexMasses({0,2,0});tool.setTensionEdges({0,1});tool.setEdgeForceMagnitudes({std::sqrt(2.0),std::sqrt(2.0)});r=tool.createForcePolyhedra();check(r.equilibrium&&tool.axialForces()[0]>0,"Tension sign and load equilibrium");
    // Deliberately unequal graph slope: AP establishes a force equilibrium;
    // reciprocation aligns the graph and both member-face normals afterwards.
    for(double weight:{0.0,.5,1.0}){
        auto initial=makeGraph({{-1,0,0},{.3f,0,.6f},{1,0,0}},{0,1,1,2});
        tool.setInputGraph(initial);tool.setSupports({0,2});tool.setVertexMasses({0,2,0});
        zForcePolyhedraSettings settings;settings.formWeight=weight;tool.setSettings(settings);
        auto before=tool.createForcePolyhedra();auto p0=points(tool);auto mesh0=tool.forceCell(1);zPointArray c0;zFnMesh(mesh0).getVertexPositions(c0);
        r=tool.solveEquilibrium(4000);print(r);auto p=points(tool);
        check(p[0].x==-1&&p[0].z==0&&p[2].x==1&&p[2].z==0,"Reciprocation preserves supports");
        if(weight==0){check(p[1].x==p0[1].x&&p[1].z==p0[1].z,"Weight zero freezes graph exactly");check(r.equilibrium,"Fixed-form compression arch resolves face areas");}
        if(weight==1){auto mesh=tool.forceCell(1);zPointArray c;zFnMesh(mesh).getVertexPositions(c);check(c.size()==c0.size(),"Weight one freezes force topology");for(size_t i=0;i<c.size();++i)check(c[i].x==c0[i].x&&c[i].y==c0[i].y&&c[i].z==c0[i].z,"Weight one freezes force positions exactly");}
        if(weight==.5){check(r.equilibrium,"Weighted 3DGS reciprocation converges");check(r.maxAngleDeviation<settings.angleTolerance,"Both face normals parallel to form edges");check(r.maxResidual<settings.equilibriumTolerance,"Physical force equilibrium");}
        inspect(tool,3);tool.reset();auto reset=points(tool);check(reset[1].x==p0[1].x&&reset[1].z==p0[1].z,"Reset restores input");
    }
    // A skew two-member node cannot balance a vertical load with fixed edges.
    auto impossible=makeGraph({{-1,0,0},{.2f,.1f,.7f},{1,0,0}},{0,1,1,2});
    tool.setInputGraph(impossible);tool.setSupports({0,2});tool.setVertexMasses({0,2,0});tool.setFormWeight(0);
    auto fixedBefore=points(tool);r=tool.solveEquilibrium(20);
    check(!r.equilibrium&&r.status=="reciprocation_incomplete","Incompatible fixed form remains an iterative solve");
    check(r.maxResidual>.01,"Incompatible form retains physical imbalance");check(points(tool)[1].y==fixedBefore[1].y,"Incompatible form stays fixed");
    auto hanging=makeGraph({{-1,0,0},{0,0,-1},{1,0,0}},{0,1,1,2});
    tool.setInputGraph(hanging);tool.setSupports({0,2});tool.setVertexMasses({0,2,0});tool.setFormWeight(0);
    r=tool.solveEquilibrium(20);check(!r.equilibrium,"Compression signs cannot silently flip to tension");
    tool.setTensionEdges({0,1});r=tool.solveEquilibrium(2000);check(r.equilibrium,"Assigned tension signs permit a hanging fixed form");
    // Force-fixed solve after editing only the form: the legacy Polytopal case.
    tool.setInputGraph(arch);tool.setSupports({0,2});tool.setVertexMasses({0,2,0});tool.setEdgeForceMagnitudes({std::sqrt(2.0),std::sqrt(2.0)});
    tool.createForcePolyhedra();auto edited=points(tool);edited[1]=zPoint(.2f,0,.7f);tool.setFormPositions(edited);tool.setFormWeight(1);
    r=tool.solveEquilibrium(2000);print(r);check(r.equilibrium,"Force-fixed reciprocation recovers the form");
    // Loaded chain, as in the user's first sketch: forces must redistribute.
    auto chain=makeGraph({{-2,0,0},{-1,0,.5f},{0,0,.8f},{1,0,.5f},{2,0,0}},{0,1,1,2,2,3,3,4});
    tool.setInputGraph(chain);tool.setSupports({0,4});tool.setVertexMasses({0,1,1,1,0});
    r=tool.solveEquilibrium(6000);print(r);check(r.equilibrium,"Loaded multi-node chain 3DGS solve");inspect(tool,5);
    auto branch=makeGraph({{1,0,0},{-.5f,.8660254f,0},{-.5f,-.8660254f,0},{0,0,1}},{0,3,1,3,2,3});
    tool.setInputGraph(branch);tool.setSupports({0,1,2});tool.setVertexMasses({0,0,0,3});tool.setEdgeForceMagnitudes({std::sqrt(2.0),std::sqrt(2.0),std::sqrt(2.0)});
    r=tool.createForcePolyhedra();print(r);check(r.equilibrium,"Spatial branch equilibrium");inspect(tool,4);check(tool.cellReport(3).forceRank==3,"Spatial tetrahedral cell");
    auto skew=makeGraph({{1,0,0},{-.5f,.8660254f,0},{-.5f,-.8660254f,0},{.1f,.15f,.7f}},{0,3,1,3,2,3});
    tool.setInputGraph(skew);tool.setSupports({0,1,2});tool.setVertexMasses({0,0,0,3});r=tool.solveEquilibrium(6000);print(r);
    check(r.equilibrium,"Skew spatial branching 3DGS solve");
    // Six balanced normals give a cube; collinear support cells receive four virtual faces.
    auto star=makeGraph({{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}},{0,1,0,2,0,3,0,4,0,5,0,6});
    tool.setInputGraph(star);tool.setSupports({1,2,3,4,5,6});tool.setEdgeForceMagnitudes({1,1,1,1,1,1});
    r=tool.createForcePolyhedra();check(r.equilibrium,"Cube force cell");inspect(tool,7);check(tool.cellFaces(0).size()==6,"Cube has six physical faces");
    zPointArray ip={{-.04f,-.06f,.02f},{1,0,0},{0,1,0},{0,0,1},{-1,-.3f,-.2f},{-.2f,-1,-.7f}};
    auto irregular=makeGraph(ip,{0,1,0,2,0,3,0,4,0,5});tool.setInputGraph(irregular);tool.setSupports({1,2,3,4,5});
    zDoubleArray target;for(int i=1;i<6;++i){auto d=ip[i]-ip[0];target.push_back(d.length());}tool.setEdgeForceMagnitudes(target);
    r=tool.createForcePolyhedra();print(r);check(r.cellsConverged,"Irregular EGI/AP converges");inspect(tool,6);int zero=0;for(auto f:tool.cellFaces(0))if(f.type==zForcePolyhedraFaceType::Zero)++zero;check(zero>0,"EGI cross-adjacencies add zero faces");
    // Corresponding-area acceptance is absolute and independent of cell geometry.
    tool.setInputGraph(chain);tool.setSupports({0,4});tool.setVertexMasses({0,1,1,1,0});
    zForcePolyhedraSettings pairSettings;pairSettings.maxIterations=1;pairSettings.areaIterationsPerStep=1;
    tool.setSettings(pairSettings);tool.createForcePolyhedra();r=tool.stepEquilibrium(1);
    const double difference=r.maxCorrespondingAreaDifference;
    check(difference>1e-10,"Nonzero paired-area difference for tolerance regression");
    pairSettings.correspondingAreaTolerance=difference*2;tool.setSettings(pairSettings);
    check(tool.report().areaPursuitResolved,"Loose corresponding-area threshold passes");
    check(!tool.report().equilibrium,"Area match alone cannot imply equilibrium");
    pairSettings.correspondingAreaTolerance=difference*.5;tool.setSettings(pairSettings);
    check(!tool.report().areaPursuitResolved,"Strict corresponding-area threshold fails");
    check(tool.report().maxCorrespondingAreaDifference==difference,"Tolerance edits preserve measured areas");
    pairSettings.correspondingAreaTolerance=-1;rejects([&]{tool.setSettings(pairSettings);});
    // No fake support reactions or success for globally unsupported loads.
    tool.setInputGraph(arch);tool.setVertexMass(1);r=tool.solveEquilibrium(20);check(!r.equilibrium,"Unsupported loaded graph cannot claim equilibrium");
    rejects([&]{tool.setSupports({-1});});rejects([&]{tool.setSupports({0,0});});rejects([&]{tool.setVertexMass(-1);});
    rejects([&]{tool.setEdgeForceMagnitudes({0,1});});rejects([&]{tool.setTensionEdges({99});});rejects([&]{tool.setVertexLoads({});});
    zForcePolyhedraSettings invalid;invalid.formWeight=1.1;rejects([&]{tool.setSettings(invalid);});rejects([&]{tool.stepEquilibrium(0);});
    auto empty=makeGraph({{0,0,0}},{});tool.setInputGraph(empty);r=tool.createForcePolyhedra();check(r.equilibrium,"Unloaded isolated node");check(tool.cellReport(0).status=="no_forces","No-force empty cell");
    std::cout<<"Force polyhedra: weighted 3DGS, virtual faces, supports, self weight, cell exports and validation passed.\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
