#include <zspace/toolsets.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <nlohmann/json.hpp>

using namespace zSpace;
void require(bool value,const char* message) { if(!value) throw std::runtime_error(message); }
zObjectGraph graph(zPointArray p,zIntArray e) { zObjectGraph g; zFnGraph(g).create(p,e); return g; }
double length(zPoint a,zPoint b) { double x=a.x-b.x,y=a.y-b.y,z=a.z-b.z; return std::sqrt(x*x+y*y+z*z); }
void checkUnroll(zObjectMesh& source,zObjectMesh& flat)
{
    zPointArray a,b; zIntArray e;
    zFnMesh(source).getVertexPositions(a); zFnMesh(flat).getVertexPositions(b); zFnMesh(source).getEdgeData(e);
    require(a.size()==b.size(),"Unroll vertex count");
    for(auto p:b) require(std::abs(p.z)<1e-6 && std::isfinite(p.x) && std::isfinite(p.y),"Flat finite coordinates");
    for(size_t i=0;i<e.size();i+=2) require(std::abs(length(a[e[i]],a[e[i+1]])-length(b[e[i]],b[e[i+1]]))<2e-4,"Unroll edge metric");
}
int main(int argc,char** argv)
{
 std::cout << std::unitbuf;
 try {
    // Deliberately shuffled connectivity: graph vertex storage is not traversal order.
    auto curve=graph({{0,0,0},{0,0,4},{0,0,2}},{2,1,0,2});
    auto profile=graph({{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0}},{2,3,0,1,3,0,1,2});
    zTsCTubes tube; tube.setDirectrix(curve); tube.setGeneratrix(profile);
    zCTubeSettings settings; settings.samples=9; tube.setSettings(settings);
    for(auto mode:{zCTubeMode::Generic,zCTubeMode::ParallelTransport})
    {
        tube.computeTube(mode); zFnMesh fn(tube.tubeMesh());
        require(fn.numVertices()==36 && fn.numPolygons()==32,"Tube topology");
        zPointArray p; fn.getVertexPositions(p);
        for(int j=0;j<9;++j) for(int k=0;k<4;++k)
        { require(std::abs(p[j*4+k].z-j*0.5)<1e-5,"Straight spline stations"); require(std::abs(std::hypot(p[j*4+k].x,p[j*4+k].y)-std::sqrt(2.0))<1e-5,"Profile radius"); }
        tube.unrollStrips(); require(tube.stripMeshes().size()==4,"Strip count");
        for(int k=0;k<4;++k) checkUnroll(tube.stripMeshes()[k],tube.unrolledMeshes()[k]);
    }
    auto bent=graph({{0,0,0},{0,0,2},{1,0,4},{2,1,6}},{0,1,1,2,2,3});
    tube.setDirectrix(bent); settings.spline=false; tube.setSettings(settings);
    for(auto mode:{zCTubeMode::Generic,zCTubeMode::ParallelTransport})
    { tube.computeTube(mode); tube.unrollStrips(); for(int k=0;k<4;++k) checkUnroll(tube.stripMeshes()[k],tube.unrolledMeshes()[k]); }
    auto loop=graph({{0,0,0},{3,0,0},{3,3,0},{0,3,0}},{0,1,1,2,2,3,3,0});
    tube.setDirectrix(loop); settings.spline=true; settings.samples=25; settings.scale=0.1; tube.setSettings(settings);
    tube.computeTube(); require(zFnMesh(tube.tubeMesh()).numVertices()==100,"Closed spline sampling");
    zPointArray sampled; zFnGraph(tube.sampledDirectrix()).getVertexPositions(sampled);
    require(length(sampled[0],sampled[24])<1e-6 && length(sampled[6],zPoint(3,0,0))<1e-6,"Periodic spline interpolation");
    // A rotated asymmetric section must retain its supplied world orientation.
    auto straight=graph({{0,0,0},{0,0,3}},{0,1});
    auto oriented=graph({{2,3,0},{3,4,0},{2,5,0},{1,4,0}},{0,1,1,2,2,3,3,0});
    tube.setDirectrix(straight); tube.setGeneratrix(oriented);
    settings.spline=false; settings.scale=1; settings.preserveOrientation=true; tube.setSettings(settings);
    for(auto mode:{zCTubeMode::Generic,zCTubeMode::ParallelTransport}) {
        tube.computeTube(mode); zPointArray q; zFnMesh(tube.tubeMesh()).getVertexPositions(q);
        require(length(q[0],zPoint(0,-1,0))<1e-6 && length(q[1],zPoint(1,0,0))<1e-6,"Preserved profile rotation");
        require(tube.frames().size()==2,"Frame count");
        require(length(tube.frames()[0].x,zPoint(float(std::sqrt(.5)),float(std::sqrt(.5)),0))<1e-5,"Frame orientation");
        tube.unrollStrips(); require(tube.stripAreas().size()==4,"Area diagnostics count");
        for(const auto& a:tube.stripAreas()) require(a.sourceArea>0 && a.deviationPercent<0.001,"Unroll area conservation");
    }
    // Optimise the apex function to pass through a larger end section of ONE tube.
    auto fitPath=graph({{0,0,0},{0,0,2},{0,0,4}},{0,1,1,2});
    zObjectGraphArray targets;
    targets.push_back(graph({{-1,-1,0},{1,-1,0},{1,1,0},{-1,1,0}},{0,1,1,2,2,3,3,0}));
    targets.push_back(graph({{-2,-2,4},{2,-2,4},{2,2,4},{-2,2,4}},{0,1,1,2,2,3,3,0}));
    tube.setDirectrix(fitPath);tube.setGeneratrices(targets);
    settings.spline=false;settings.fitPlaneNormals=false;settings.checkFitGradient=true;tube.setSettings(settings);tube.computeTube();
    require(tube.fitReport().gradientCheckError<1e-5,"Apex analytic gradient matches central differences");
    require(tube.fitReport().active && tube.fitReport().converged,"Target-section fit convergence");
    require(tube.fitReport().rms<tube.fitReport().initialRms*0.01,"Target fit improves error");
    require(tube.fitReport().apexScales.size()==2 && tube.fitReport().apexScales[0]>1.1,"Per-segment apex values optimised");
    require(zFnMesh(tube.tubeMesh()).numPolygons()==8,"Single fitted tube topology");
    tube.unrollStrips();for(const auto& a:tube.stripAreas())require(a.deviationPercent<.001,"Fitted tube unroll area");
    // A picked vertex becomes the geometric anchor in both construction modes.
    tube.setDirectrix(straight);tube.setGeneratrix(oriented);tube.setGeneratrixOrigin(0,zPoint(3,4,0));
    settings.fitPlaneNormals=false;tube.setSettings(settings);
    for(auto mode:{zCTubeMode::Generic,zCTubeMode::ParallelTransport}) {
        tube.computeTube(mode);zPointArray points;zFnMesh(tube.tubeMesh()).getVertexPositions(points);
        require(length(points[1],zPoint(0,0,0))<1e-6,"Picked vertex lies on directrix");
        require(length(points[0],zPoint(-1,-1,0))<1e-6,"Origin change preserves section shape");
    }
    auto section=[](float z,float r,float tilt=0.f) {
        float c=std::cos(tilt),s=std::sin(tilt);
        return graph({{-r*c,-r,z+r*s},{r*c,-r,z-r*s},{r*c,r,z-r*s},{-r*c,r,z+r*s}},{0,1,1,2,2,3,3,0});
    };
    auto path4=graph({{0,0,0},{0,0,1},{0,0,2},{0,0,3},{0,0,4}},{0,1,1,2,2,3,3,4});
    zObjectGraphArray varying={section(0,1),section(1,2),section(2,.7f),section(4,2)};
    tube.setDirectrix(path4);tube.setGeneratrices(varying);settings.fitIterations=100;
    settings.apexSmoothnessWeight=0;tube.setSettings(settings);tube.computeTube();auto exact=tube.fitReport();
    auto roughness=[](const zDoubleArray& a){double e=0;for(size_t j=1;j+1<a.size();++j){double d=a[j-1]-2*a[j]+a[j+1];e+=d*d;}return e;};
    settings.apexSmoothnessWeight=10;tube.setSettings(settings);tube.computeTube();auto smooth=tube.fitReport();
    require(roughness(smooth.apexScales)<roughness(exact.apexScales)*.5,"Apex smoothness reduces second differences");
    require(smooth.objective<smooth.initialObjective,"Combined weighted objective improves");
    settings.sectionWeights={1,1000,1,1};tube.setSettings(settings);tube.computeTube();
    require(tube.fitReport().sectionRms[1]<smooth.sectionRms[1],"Per-target section weight prioritises target");
    settings.sectionWeights.clear();settings.apexSmoothnessWeight=0;settings.targetSectionWeight=0;
    settings.fitPlaneNormals=true;settings.planeSmoothnessWeight=10;settings.targetPlaneWeight=0;
    zObjectGraphArray tilted={section(0,1),section(4,1,.5f)};tube.setGeneratrices(tilted);tube.setSettings(settings);tube.computeTube();
    require(tube.fitReport().gradientCheckError<1e-5,"Normal and smoothness gradients match central differences");
    double relaxedAngle=tube.fitReport().planeAngleDegrees[1];
    settings.targetPlaneWeight=1000;tube.setSettings(settings);tube.computeTube();
    require(tube.fitReport().planeAngleDegrees[1]<relaxedAngle*.1,"Plane weight prioritises target orientation");
    require(tube.fitReport().objective<tube.fitReport().initialObjective,"Plane objective improves");
    // All objective terms on an asymmetric curved path, including backward propagation
    // from an interior seed and multiple independently weighted target sections.
    auto gradientPath=graph({{0,0,-1},{0,0,0},{.2f,0,1},{.4f,.1f,2},{.1f,0,3},{0,0,4}},{0,1,1,2,2,3,3,4,4,5});
    zObjectGraphArray gradientTargets={section(0,1),section(2,1.4f),section(4,.8f,.15f)};
    tube.setDirectrix(gradientPath);tube.setGeneratrices(gradientTargets);
    settings.targetSectionWeight=1;settings.targetPlaneWeight=.7;settings.apexSmoothnessWeight=.03;settings.planeSmoothnessWeight=.02;
    settings.sectionWeights={1,.5,2};settings.checkFitGradient=true;tube.setSettings(settings);tube.computeTube();
    std::cout<<"Analytic gradient relative error: "<<tube.fitReport().gradientCheckError<<'\n';
    require(tube.fitReport().gradientCheckError<1e-5,"Full weighted objective analytic gradient");
    require(tube.fitReport().objective<tube.fitReport().initialObjective,"L-BFGS-B reduces combined objective");
    settings.targetSectionWeight=-1;bool badWeight=false;try{tube.setSettings(settings);}catch(const std::invalid_argument&){badWeight=true;}
    require(badWeight,"Negative objective weights rejected");
    auto nonplanarProfile=graph({{300,100,0},{302,100,0},{302,102,.1f},{300,102,0}},{0,1,1,2,2,3,3,0});
    bool warpedRejected=false;try{tube.setGeneratrix(nonplanarProfile);}catch(const std::invalid_argument&){warpedRejected=true;}
    require(warpedRejected,"World-coordinate tolerance still rejects nonplanar profiles");
    auto branch=graph({{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{0,1,0,2,0,3});
    bool rejected=false; try { tube.setDirectrix(branch); } catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"Reject branched graph");
    auto open=graph({{-1,0,0},{1,0,0}},{0,1});
    settings=zCTubeSettings{};settings.spline=false;settings.preserveOrientation=true;tube.setSettings(settings);
    tube.setDirectrix(straight);tube.setGeneratrix(open);
    for(auto mode:{zCTubeMode::Generic,zCTubeMode::ParallelTransport}) {
        tube.computeTube(mode);require(zFnMesh(tube.tubeMesh()).numPolygons()==1,"Open line produces one strip without seam face");
        tube.unrollStrips();require(tube.stripMeshes().size()==1,"Two-point generatrix has one strip");
        require(tube.stripAreas()[0].deviationPercent<.001,"Open strip unroll conserves area");
    }
    zObjectGraphArray openTargets={open,graph({{-2,0,4},{2,0,4}},{0,1})};
    tube.setDirectrix(fitPath);tube.setGeneratrices(openTargets);settings.fitPlaneNormals=false;tube.setSettings(settings);tube.computeTube();
    require(tube.fitReport().converged,"Open line target fit converges");require(zFnMesh(tube.tubeMesh()).numPolygons()==2,"Open fit has no closing face");
    tube.unrollStrips();require(tube.stripMeshes().size()==1,"Open fit unrolls one strip");
    zObjectGraphArray openChains={graph({{-1,0,0},{0,1,0},{1,0,0}},{0,1,1,2}),graph({{-2,0,4},{0,2,4},{2,0,4}},{0,1,1,2})};
    tube.setGeneratrices(openChains);tube.computeTube();require(tube.fitReport().converged,"Open polyline target fit converges");
    require(zFnMesh(tube.tubeMesh()).numPolygons()==4,"Three-point open profile produces two strips");
    tube.unrollStrips();require(tube.stripMeshes().size()==2,"Open chain strip count");
    for(size_t k=0;k<2;++k)checkUnroll(tube.stripMeshes()[k],tube.unrolledMeshes()[k]);
    auto tSection=graph({{0,0,0},{-1,0,0},{1,0,0},{0,-1,0}},{0,1,0,2,0,3});
    tube.setDirectrix(fitPath);tube.setGeneratrix(tSection);tube.setGeneratrixOrigin(0,{0,0,0});
    for(auto mode:{zCTubeMode::Generic,zCTubeMode::ParallelTransport}){
        tube.computeTube(mode);require(zFnMesh(tube.tubeMesh()).numPolygons()==6,"T section sweeps its three graph edges");
        tube.unrollStrips();require(tube.stripMeshes().size()==2,"T crossbar and stem produce two strips");
        for(size_t k=0;k<2;++k){checkUnroll(tube.stripMeshes()[k],tube.unrolledMeshes()[k]);require(tube.stripAreas()[k].deviationPercent<.001,"T strip area conserved");}
    }
    // Offset sections need no physical contact with the directrix. Vertex IDs and line directions may differ.
    zObjectGraphArray tTargets={graph({{3,2,0},{2,2,0},{4,2,0},{3,1,0}},{0,1,0,2,0,3}),graph({{6,2,4},{8,4,4},{6,4,4},{4,4,4}},{2,3,1,2,0,2})};
    tube.setGeneratrices(tTargets);tube.setSettings(settings);tube.computeTube();
    require(tube.fitReport().converged,"Offset T targets with permuted vertex IDs fit");
    tube.unrollStrips();require(tube.stripMeshes().size()==2,"Fitted T combines collinear flange edges");
    for(size_t k=0;k<2;++k)checkUnroll(tube.stripMeshes()[k],tube.unrolledMeshes()[k]);
    std::reverse(tTargets.begin(),tTargets.end());tube.setGeneratrices(tTargets);tube.computeTube();
    require(tube.fitReport().converged,"T target topology survives reversed section order");
    auto otherBranch=graph({{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0}},{0,1,0,2,0,3,0,4});
    auto subdividedT=graph({{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,2,0}},{0,1,0,2,0,3,3,4});
    zObjectGraphArray mismatched={otherBranch,subdividedT};rejected=false;try{tube.setGeneratrices(mismatched);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"Different branch connectivity rejected");
    rejected=false;try{tube.setGeneratrix(branch);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"Nonplanar section graph rejected");
    auto splitLine=graph({{-1,0,0},{0,0,0},{1,0,0}},{1,0,2,1});tube.setGeneratrix(splitLine);tube.computeTube();tube.unrollStrips();
    require(tube.stripMeshes().size()==1,"Collinear chain edges unroll together");checkUnroll(tube.stripMeshes()[0],tube.unrolledMeshes()[0]);
    auto ySection=graph({{0,0,0},{-1,0,0},{.5f,1,0},{.5f,-1,0}},{0,1,0,2,0,3});tube.setGeneratrix(ySection);tube.computeTube();tube.unrollStrips();
    require(tube.stripMeshes().size()==3,"Noncollinear branches keep separate strips");
    auto longPath=graph({{0,0,0},{0,0,50},{0,0,100}},{0,1,1,2});
    zObjectGraphArray nearEnd={graph({{0,0,.0005f},{-1,0,.0005f},{1,0,.0005f},{0,-1,.0005f}},{0,1,0,2,0,3}),graph({{0,0,100},{-2,0,100},{2,0,100},{0,-2,100}},{0,1,0,2,0,3})};
    tube.setDirectrix(longPath);tube.setGeneratrices(nearEnd);settings.apexSmoothnessWeight=.01;settings.planeSmoothnessWeight=.01;tube.setSettings(settings);tube.computeTube();
    require(tube.frames().size()==3,"Near-endpoint target replaces display station instead of creating a sliver");
    require(std::abs(tube.frames()[0].origin.z-.0005)<1e-7,"Exact target plane station retained");
    tube.unrollStrips();require(tube.stripMeshes().size()==2,"Near-endpoint T unrolls both strips");
    // Inward flange tips follow a variable central seam while the crossbar stays straight.
    auto fixedPath=graph({{0,0,0},{0,0,1},{0,0,2}},{0,1,1,2});
    auto fixedT=graph({{0,0,0},{-.2f,0,0},{.2f,0,0},{0,-.4f,0}},{1,0,0,2,0,3});
    std::vector<zCTubeFrame> fixedFrames={{{0,0,0},{1,0,0},{0,1,0},{0,0,1}},{{0,0,1},{1,0,0},{0,1,0},{0,0,1}},{{0,0,2},{1,0,0},{0,1,0},{0,0,1}}};
    zPointArray seam={{.3f,0,0},{.6f,0,1},{.9f,0,2}};
    tube.setDirectrix(fixedPath);tube.setGeneratrix(fixedT);settings.spline=false;tube.setSettings(settings);
    tube.computeTubeWithFrames(fixedFrames,seam,2);
    zPointArray closedPositions;zFnMesh(tube.tubeMesh()).getVertexPositions(closedPositions);
    for(int j=0;j<3;++j)require(length(closedPositions[6+2*j+1],seam[j])<1e-7,"Inward flange reaches exact central seam");
    tube.computeStrips();tube.unrollStrips();require(tube.stripMeshes().size()==2,"Variable inward flange stays one flange and one stem strip");
    // Nonplanar strip exercises the core triangulation fallback.
    zObjectMesh warped,flat; zPointArray p={{0,0,0},{2,0,0},{2,1,1},{0,1,0}}; zIntArray c={4}, f={0,1,2,3};
    zFnMesh(warped).create(p,c,f); zFnMesh(warped).unroll(flat); checkUnroll(warped,flat);
    rejected=false; try { zFnMesh(warped).unroll(flat,false); } catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"Reject nonplanar face when triangulation disabled");
    zObjectMesh tetra; p={{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    c={3,3,3,3}; f={0,2,1,0,1,3,1,2,3,2,0,3}; zFnMesh(tetra).create(p,c,f);
    rejected=false; try { zFnMesh(tetra).unroll(flat); } catch(const std::invalid_argument&) { rejected=true; }
    require(rejected,"Closed nondevelopable mesh needs seams");
    // Coupled seams may move while all boundary anchor positions stay fixed.
    {
        std::vector<zCTubeBranchInput> branches(3);
        for(int b=0;b<3;++b){double angle=b*2*3.141592653589793/3;zPoint x(float(-std::cos(angle)),float(-std::sin(angle)),0),y(float(std::sin(angle)),float(-std::cos(angle)),0);for(int j=0;j<5;++j){double radius=1+.08*(b+1)*std::sin(j*.6);branches[b].frames.push_back({zPoint(float(radius*std::cos(angle)),float(radius*std::sin(angle)),float(j)),x,y,zPoint(0,0,1)});branches[b].centralTargets.push_back(zPoint(float(.04*std::sin(j*.7)),float(.03*std::cos(j*.5)),float(j)));}}
        zTsCTubeBranching optimiser;zCTubeBranchSettings options;options.fixedDirectrices=true;options.iterations=100;options.tolerance=1e-5;optimiser.optimise(branches,options);
        require(optimiser.report().gradientCheckError<1e-5,"Coupled seam gradient agrees with finite differences");require(optimiser.report().converged,"Fixed boundary seam optimisation converges");
        for(size_t b=0;b<3;++b){const auto& result=optimiser.branches()[b];zPointArray origins,seed;zIntArray edges;for(size_t j=0;j<result.frames.size();++j){require(length(result.frames[j].origin,branches[b].frames[j].origin)==0,"Central optimisation preserves boundary anchors exactly");origins.push_back(result.frames[j].origin);if(j)edges.insert(edges.end(),{int(j-1),int(j)});}seed.assign(result.sections.begin(),result.sections.begin()+4);auto path=graph(origins,edges),section=graph(seed,{1,0,0,2,0,3});tube.setDirectrix(path);tube.setGeneratrix(section);tube.computeTubeWithSections(result.frames,result.sections);tube.unrollStrips();require(tube.stripMeshes().size()==2,"Fixed boundary optimisation retains two strips");for(size_t i=0;i<2;++i)checkUnroll(tube.stripMeshes()[i],tube.unrolledMeshes()[i]);}
        std::cout<<"Fixed boundary central optimisation: seam gap "<<optimiser.report().seamError<<", gradient error "<<optimiser.report().gradientCheckError<<"\n";
    }
    if(argc>1) {
        std::ifstream stream(argv[1]);require(bool(stream),"Open closed-flange fixture");nlohmann::json input;stream>>input;
        auto positions=[](const nlohmann::json& a){zPointArray p;for(size_t i=0;i<a.size();i+=3)p.emplace_back(a[i].get<float>(),a[i+1].get<float>(),a[i+2].get<float>());return p;};
        auto xyz=[](const nlohmann::json& a){return zPoint(a[0].get<float>(),a[1].get<float>(),a[2].get<float>());};
        for(const auto& branch:input) {
            auto path=graph(positions(branch["directrix"]["vertexPositions"]),branch["directrix"]["edgeConnects"].get<zIntArray>());
            auto section=graph(positions(branch["generatrix"]["vertexPositions"]),branch["generatrix"]["edgeConnects"].get<zIntArray>());
            std::vector<zCTubeFrame> frames;for(const auto& f:branch["frames"])frames.push_back({xyz(f["origin"]),xyz(f["x"]),xyz(f["y"]),xyz(f["normal"])});
            zPointArray targets;for(const auto& p:branch["flangeTargets"])targets.push_back(xyz(p));
            tube.setDirectrix(path);tube.setGeneratrix(section);tube.computeTubeWithFrames(frames,targets,2);
            zPointArray p;zFnMesh(tube.tubeMesh()).getVertexPositions(p);
            for(size_t j=0;j<frames.size();++j)require(length(p[2*frames.size()+2*j+1],targets[j])<1e-7,"Actual Y mesh tips meet central targets");
            tube.computeStrips();require(tube.stripMeshes().size()==2,"Actual closed Y flanges retain two strips");
            try {tube.unrollStrips();std::cout<<"Closed Y branch: "<<frames.size()<<" stations; exact unroll accepted\n";}
            catch(const std::invalid_argument& error) {
                require(std::string(error.what()).find("without a seam or metric distortion")!=std::string::npos,"Only nondevelopability may reject closed flange unroll");
                require(tube.unrolledMeshes().empty(),"Rejected unroll must not expose a partially flattened result");
                std::cout<<"Closed Y branch: "<<frames.size()<<" stations; needs a cut or distortion\n";
            }
        }
        for(size_t first=0;first<input.size();first+=3) {
            std::vector<zCTubeBranchInput> branches;
            for(size_t b=first;b<first+3;++b){branches.emplace_back();for(const auto& f:input[b]["frames"])branches.back().frames.push_back({xyz(f["origin"]),xyz(f["x"]),xyz(f["y"]),xyz(f["normal"])});for(const auto& p:input[b]["flangeTargets"])branches.back().centralTargets.push_back(xyz(p));}
            zTsCTubeBranching optimiser;zCTubeBranchSettings options;options.iterations=50;optimiser.optimise(branches,options);
            require(optimiser.report().gradientCheckError<1e-5,"Branch analytic gradient agrees with finite differences");
            for(size_t b=0;b<3;++b){const auto& result=optimiser.branches()[b];zPointArray origins,seed;zIntArray edges;for(size_t j=0;j<result.frames.size();++j){origins.push_back(result.frames[j].origin);if(j)edges.insert(edges.end(),{int(j-1),int(j)});require(length(result.sections[4*j+2],branches[b].centralTargets[j])<1e-7,"Optimisation preserves exact shared seam");}seed.assign(result.sections.begin(),result.sections.begin()+4);auto path=graph(origins,edges),section=graph(seed,{1,0,0,2,0,3});tube.setDirectrix(path);tube.setGeneratrix(section);tube.computeTubeWithSections(result.frames,result.sections);tube.computeStrips();tube.unrollStrips();require(tube.stripMeshes().size()==2,"Optimised joined tube has two exactly unrolled strips");for(size_t i=0;i<2;++i)checkUnroll(tube.stripMeshes()[i],tube.unrolledMeshes()[i]);}
            std::cout<<"Branch optimiser: RMS "<<optimiser.report().rmsDisplacement<<", max "<<optimiser.report().maxDisplacement<<", gradient error "<<optimiser.report().gradientCheckError<<"; all six strips unroll exactly\n";
        }
    }
    std::cout<<"C-Tubes generation and unroll checks passed\n";
    return 0;
 } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
