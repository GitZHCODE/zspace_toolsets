#include "../../src/zToolsets/fabrication/zTsStairsLayout.h"
#include <zspace/toolsets.h>
#include <iostream>
#include <cmath>
#include <stdexcept>
using namespace zSpace;
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
template<class F>void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"Expected rejection");}
int main(){try{
    zObjectGraph g;zPointArray p={{-5,0,0},{0,0,2},{5,0,0},{0,5,0}};zIntArray e={0,1,1,2,1,3};zFnGraph(g).create(p,e);
    zTsParametricStairs t;zParametricStairsSettings s;s.width=1;s.divisions=4;s.vertexWidths={1,2,1,1};t.update(g,s);
    check(t.supports().size()==9,"All three leaf caps supported");check(t.centerGraphMeshVertices().size()==4,"Graph centers mapped");
    zPointArray before;zFnMesh(t.centerMesh()).getVertexPositions(before);
    check(before.size()==40,"Connected T ribbon vertex count");check(zFnMesh(t.centerMesh()).numPolygons()==24,"Two quad strips per arm");
    for(int v:t.supports())check(t.meshVertexGraphVertices()[v]!=1,"Junction is free");
    for(int i=0;i<4;++i){auto q=before[t.centerGraphMeshVertices()[i]];check(q.x==p[i].x&&q.y==p[i].y&&q.z==p[i].z,"Centers match source");}
    // The terminal cap is one model unit wide; at the straight-through junction it is two.
    double minY=1e9,maxY=-1e9;for(size_t i=0;i<before.size();++i)if(t.meshVertexGraphVertices()[i]==0){minY=std::min(minY,double(before[i].y));maxY=std::max(maxY,double(before[i].y));}
    check(std::abs(maxY-minY-1)<1e-6,"Per-vertex width respected");
    auto initial=t.report();t.stepBestFit(5);zPointArray after;zFnMesh(t.resultMesh()).getVertexPositions(after);
    for(size_t i=0;i<before.size();++i)check(after[i].x==before[i].x&&after[i].y==before[i].y,"XY fixed by default");
    for(int v:t.supports())check(after[v].z==before[v].z,"Support elevation fixed");
    for(double q:t.forceDensities())check(std::isfinite(q)&&q<=0,"Compression-only force densities");
    check(std::isfinite(initial.residualNorm)&&std::isfinite(t.report().residualNorm),"Finite solver diagnostics");
    auto invalid=s;invalid.width=-1;rejects([&]{t.update(g,invalid);});check(t.supports().size()==9,"Invalid update preserves result");
    invalid=s;invalid.vertexWidths={1};rejects([&]{t.update(g,invalid);});invalid=s;invalid.divisions=1;rejects([&]{t.update(g,invalid);});
    zObjectGraph loop;zPointArray lp={{0,0,0},{2,0,0},{1,2,0}};zIntArray le={0,1,1,2,2,0};zFnGraph(loop).create(lp,le);rejects([&]{t.update(loop);});
    zObjectGraph line;zPointArray l={{0,0,0},{4,0,1},{8,0,0}};zIntArray lc={0,1,1,2};auto lineSettings=s;lineSettings.vertexWidths.clear();zFnGraph(line).create(l,lc);t.update(line,lineSettings);check(t.supports().size()==6,"Straight chain joins");
    l[1].z=2;zFnGraph(line).create(l,lc);t.update(line,lineSettings);zFnMesh(t.centerMesh()).getVertexPositions(after);check(after[t.centerGraphMeshVertices()[1]].z==2,"Live update resets reference elevation");
    // New pipeline: coarse centre graph edges, all-face division, smoothing BEFORE FDM.
    zParametricStairsSettings refined;refined.numDivs=0;refined.numSmoothDivs=0;t.update(g,refined);
    check(zFnMesh(t.lowPolyMesh()).numPolygons()==6&&zFnMesh(t.lowPolyMesh()).numVertices()==13,"Coarse central graph ribbon");
    zIntArray coarseEdges;zFnMesh(t.lowPolyMesh()).getEdgeData(coarseEdges);
    for(size_t i=0;i<e.size();i+=2){int a=t.centerGraphMeshVertices()[e[i]],b=t.centerGraphMeshVertices()[e[i+1]];bool found=false;for(size_t j=0;j<coarseEdges.size();j+=2)if((coarseEdges[j]==a&&coarseEdges[j+1]==b)||(coarseEdges[j]==b&&coarseEdges[j+1]==a))found=true;check(found,"Coarse mesh includes every central graph edge");}
    refined.numDivs=1;t.update(g,refined);check(zFnMesh(t.subdividedMesh()).numPolygons()==12,"Each quad splits into two along the arm");check(t.supports().size()==9,"Directional division retains original cap supports");
    check(zFnMesh(t.subdividedMesh()).numVertices()==22,"Only three vertices per new transverse station");
    zObjectMesh expected=t.subdividedMesh();zFnMesh(expected).smoothMesh(1,false);zPointArray expectedPositions;zFnMesh(expected).getVertexPositions(expectedPositions);
    refined.numSmoothDivs=1;t.update(g,refined);check(zFnMesh(t.centerMesh()).numPolygons()==48,"Smooth subdivision before FDM");check(t.supports().size()==15,"Smooth cap points inherit supports");
    zFnMesh(t.centerMesh()).getVertexPositions(before);zFnMesh(t.resultMesh()).getVertexPositions(after);
    check(before.size()==expectedPositions.size(),"Core Catmull-Clark vertex count");bool moved=false;for(size_t i=0;i<before.size();++i){check((before[i]-expectedPositions[i]).length()<1e-5,"Smoothed input matches core smoothing");check((before[i]-after[i]).length()<1e-6,"FDM starts from smoothed input");}
    zPointArray flatDivided;zFnMesh(t.subdividedMesh()).getVertexPositions(flatDivided);for(size_t i=0;i<flatDivided.size();++i)if((flatDivided[i]-before[i]).length()>1e-5)moved=true;check(moved,"Smoothing changes shape before solve");
    t.stepBestFit(3);zFnMesh(t.resultMesh()).getVertexPositions(after);for(size_t i=0;i<before.size();++i)check(before[i].x==after[i].x&&before[i].y==after[i].y,"FDM fixes smoothed XY");for(int id:t.supports())check(before[id].z==after[id].z,"Refined supports fixed");
    t.resetBestFit();zFnMesh(t.resultMesh()).getVertexPositions(after);for(size_t i=0;i<before.size();++i)check((before[i]-after[i]).length()<1e-6,"Reset restores smoothed input");
    auto invalidRefinement=refined;invalidRefinement.numDivs=-1;rejects([&]{t.update(g,invalidRefinement);});check(zFnMesh(t.centerMesh()).numPolygons()==48,"Failed refinement preserves solve");
    auto directional=refined;directional.numDivs=2;directional.numSmoothDivs=0;t.update(line,directional);zFnMesh(t.subdividedMesh()).getVertexPositions(after);
    check(zFnMesh(t.subdividedMesh()).numPolygons()==16&&after.size()==27,"Four longitudinal segments, two faces across width");
    for(auto point:after)check(point.y==-.5f||point.y==0||point.y==.5f,"No extra lengthwise lanes from face divisions");
    check(t.supports().size()==6,"End caps not divided across width");
    refined.numDivs=2;refined.numSmoothDivs=4;t.update(g,refined);
    check(zFnMesh(t.centerMesh()).numPolygons()==6144,"Requested directional 2 face / 4 smooth divisions accepted");
    check(t.report().activeEdges>1024,"Old edge limit remains removed");
    check(t.report().rank==-1&&std::isfinite(t.report().residualNorm),"Sparse initialization with unmeasured rank");
    zFnMesh(t.centerMesh()).getVertexPositions(before);auto largeReport=t.stepBestFit(1);zFnMesh(t.resultMesh()).getVertexPositions(after);
    check(std::isfinite(largeReport.residualNorm)&&largeReport.iteration>0,"Large FDM advances a finite accepted step");
    for(size_t i=0;i<before.size();++i)check(before[i].x==after[i].x&&before[i].y==after[i].y,"Large FDM preserves XY");
    for(int id:t.supports())check(before[id].z==after[id].z,"Large FDM preserves cap supports");
    for(double q:t.forceDensities())check(std::isfinite(q)&&q<=0,"Large FDM preserves compression assignments");
    t.resetBestFit();check(t.report().iteration==0,"Large FDM reset");
    // Independent equilibrium oracle for the sparse initializer: a sampled
    // paraboloid grid has uniform q=-1 and mass=4 under gravity (0,0,-1).
    const int size=33;zObjectMesh grid;zPointArray gp;zIntArray gc,gf,gs;
    for(int y=0;y<size;++y)for(int x=0;x<size;++x){float a=float(x-size/2),b=float(y-size/2);gp.emplace_back(a,b,600-a*a-b*b);if(x==0||y==0||x==size-1||y==size-1)gs.push_back(y*size+x);}
    for(int y=0;y<size-1;++y)for(int x=0;x<size-1;++x){int a=y*size+x;gc.push_back(4);gf.insert(gf.end(),{a,a+1,a+1+size,a+size});}
    zFnMesh(grid).create(gp,gc,gf);zTsVault analytic;analytic.setInputMesh(grid);analytic.setSupports(gs);analytic.setVertexMass(4);
    zVaultBestFitSettings analyticSettings;analyticSettings.fixXY=true;analytic.setBestFitSettings(analyticSettings);auto analyticReport=analytic.computeBestFitForceDensities();
    check(analyticReport.activeEdges>1024&&analyticReport.rank==-1,"Analytic mesh uses sparse path");
    std::cout<<"Sparse analytic residual "<<analyticReport.residualNorm<<" max "<<analyticReport.maxResidual<<"\n";
    check(analyticReport.residualNorm/(4*std::sqrt(double(analyticReport.freeVertices)))<1e-3&&analyticReport.maxResidual/4<1e-3,"Sparse fit matches analytic grid loads within 0.1 percent");
    check(!analyticReport.converged,"Inner iteration budget does not falsely claim outer convergence");
    // Fixed q=-0.2, load=-2 gives z=5. Attraction 0.5 gives z=3;
    // attraction 1 retains the smoothed z=1 reference.
    zObjectMesh arch;zPointArray ap={{-1,0,0},{1,0,0},{0,0,1}};zIntArray ac={3},af={0,1,2};zFnMesh(arch).create(ap,ac,af);
    zTsVault snap;snap.setInputMesh(arch);snap.setSupports({0,1});snap.setVertexMass(2);zVaultBestFitSettings ss;ss.fixXY=true;ss.qMin=.2;ss.qMax=.2;snap.setBestFitSettings(ss);snap.computeBestFitForceDensities();
    auto fixedQ=snap.forceDensities();snap.snapFixedDensities(0);zFnMesh(snap.resultMesh()).getVertexPositions(after);check(std::abs(after[2].z-5)<1e-5&&snap.report().converged,"Direct equilibrium arch");check(snap.forceDensities()==fixedQ,"Snap retains densities");
    snap.resetBestFit();snap.snapFixedDensities(.5);zFnMesh(snap.resultMesh()).getVertexPositions(after);check(std::abs(after[2].z-3)<1e-5&&!snap.report().converged,"Attraction trades equilibrium for reference proximity");
    snap.snapFixedDensities(1);zFnMesh(snap.resultMesh()).getVertexPositions(after);check(after[2].z==1,"Weight one fixes reference geometry");check(after[0].z==0&&after[1].z==0,"Snap retains supports");
    rejects([&]{snap.snapFixedDensities(-1);});check(snap.previousBestFit(),"Snap rewind available");
    ap[2]=zPoint(.3f,.4f,1);zFnMesh(arch).create(ap,ac,af);
    snap.setInputMesh(arch);snap.setSupports({0,1});snap.setVertexMass(2);ss.fixXY=false;ss.planConstraint=zVaultPlanConstraint::Free;snap.setBestFitSettings(ss);snap.computeBestFitForceDensities();
    snap.snapFixedDensities(0);zFnMesh(snap.resultMesh()).getVertexPositions(after);check(std::abs(after[2].x)<1e-6&&std::abs(after[2].y)<1e-6&&std::abs(after[2].z-5)<1e-5,"Unlocked snap solves XYZ equilibrium");
    snap.snapFixedDensities(.5);zFnMesh(snap.resultMesh()).getVertexPositions(after);check(std::abs(after[2].x-.15)<1e-6&&std::abs(after[2].y-.2)<1e-6&&std::abs(after[2].z-3)<1e-5,"Unlocked attraction uses all free coordinates");
    zParametricStairsSettings refineSettings;refineSettings.numDivs=0;refineSettings.numSmoothDivs=0;t.update(g,refineSettings);t.stepBestFit();
    zFnMesh(t.resultMesh()).getVertexPositions(before);auto pinned=t.supports();int faceCount=zFnMesh(t.resultMesh()).numPolygons();
    t.smoothResult(1);zFnMesh(t.resultMesh()).getVertexPositions(after);check(zFnMesh(t.resultMesh()).numPolygons()==4*faceCount,"Post smooth changes working topology");
    check(t.supports().size()>pinned.size(),"Post smooth adds cap supports");for(int v:pinned)check((after[v]-before[v]).length()<1e-6,"Post smooth preserves support anchors");
    check(t.forceDensities().size()==t.edgeConnects().size()/2,"Refined solver densities match current edges");t.stepBestFit();t.resetBestFit();check(zFnMesh(t.resultMesh()).numPolygons()==4*faceCount,"Subsequent solve and reset retain refined topology");
    zObjectGraph walk;zPointArray wp={{0,0,0},{4,0,.7f},{8,0,1.4f}};zIntArray we={0,1,1,2};zFnGraph(walk).create(wp,we);
    zParametricStairsSettings stairSettings;stairSettings.generateSteps=true;stairSettings.riserHeight=.2;stairSettings.treadDepth=.3;stairSettings.numSmoothDivs=0;
    t.update(walk,stairSettings);check(t.flightRiserCounts()==zIntArray{7},"Integer stair riser count");check(std::abs(t.flightActualRisers()[0]-.2)<1e-6,"Uniform actual riser");
    check(t.minimumClearance()>=stairSettings.clearance-1e-5,"Initial faces clear tread undersides");check(zFnMesh(t.centerMesh()).numPolygons()==8,"Base mesh follows numDivs topology without stair divisions");auto kinds=t.stairFaceTypes();check(std::count(kinds.begin(),kinds.end(),1)==0,"No vertical riser surfaces in stair mesh");check(std::count(kinds.begin(),kinds.end(),0)>=12&&std::count(kinds.begin(),kinds.end(),2)>0,"Open tread and landing surface panels");
    check(zFnGraph(t.walkingGraph()).numVertices()==3&&zFnGraph(t.walkingGraph()).numEdges()==2,"Editable graph keeps original nodes and edges");
    t.snapFixedDensities(0);check(t.minimumClearance()>=stairSettings.clearance-1e-5,"Snap respects stair clearance");t.stepBestFit(2);check(t.minimumClearance()>=stairSettings.clearance-1e-5,"Play respects stair clearance");
    t.smoothResult(1);check(t.minimumClearance()>=stairSettings.clearance-1e-5,"Post smoothing respects stair clearance");
    zPointArray stable;zFnMesh(t.resultMesh()).getVertexPositions(stable);stairSettings.treadDepth=3;rejects([&]{t.update(walk,stairSettings);});zFnMesh(t.resultMesh()).getVertexPositions(after);check(after.size()==stable.size()&&(after.front()-stable.front()).length()<1e-6,"Invalid flight preserves previous layout");
    // A stair footprint can intersect a structural triangle even when every vertex is outside it.
    zStairsLayoutData envelope;envelope.footprints.push_back({{{-1,1,.1f},{5,1,.1f},{5,2,.1f},{-1,2,.1f}},.1});
    zPointArray crossing={{0,0,.2f},{4,0,.2f},{0,4,.2f}};zIntArray triangleCounts={3},triangleIds={0,1,2};
    check(stairsClearance(crossing,triangleCounts,triangleIds,envelope,0)<-.09,"Clearance catches edge-crossing face intersections");
    lowerStairsBelowEnvelope(crossing,triangleCounts,triangleIds,{},envelope,0,.05);check(stairsClearance(crossing,triangleCounts,triangleIds,envelope,0)>=.05-1e-6,"Face envelope correction removes intersection");
    stairSettings.treadDepth=.3;stairSettings.landingPolygons={{{-.7f,-1,0},{.7f,-1,0},{.7f,1,0},{-.7f,1,0}},{{7.3f,-1,1.4f},{8.7f,-1,1.4f},{8.7f,1,1.4f},{7.3f,1,1.4f}}};
    t.update(walk,stairSettings);check(t.minimumClearance()>=.05-1e-5,"Rhino landing outlines form valid stair envelopes");
    zObjectGraph stacked;zPointArray sp={{0,0,0},{8,0,1.4f},{0,0,4},{8,0,5.4f}};zIntArray se={0,1,2,3};zFnGraph(stacked).create(sp,se);
    stairSettings.landingPolygons.clear();t.update(stacked,stairSettings);check(t.flightRiserCounts().size()==2,"Stacked flights retain separate ownership");check(t.minimumClearance()>=.05-1e-5,"Upper structural flight may pass above lower steps");t.stepBestFit();check(t.minimumClearance()>=.05-1e-5,"Stacked clearance remains protected during Play");
    { zObjectGraph levelGraph;zPointArray lp={{0,0,0},{2,0,.68f}};zIntArray le={0,1};zFnGraph(levelGraph).create(lp,le);
    zParametricStairsSettings os;os.generateSteps=true;os.numSmoothDivs=0;
    auto keep=zTsParametricStairs::optimiseGraph(levelGraph,os,.25,.35,.15,.18);zPointArray levels;zFnGraph(keep.graph).getVertexPositions(levels);
    check(keep.adjustedLandings==0&&levels[1].z==lp[1].z&&keep.riserCounts==zIntArray{4},"Domain fit first preserves current landing levels");
    lp[1].z=1.4f;zFnGraph(levelGraph).create(lp,le);os.landingPolygons={{{1.5f,-.5f,1.4f},{2.5f,-.5f,1.4f},{2.5f,.5f,1.4f},{1.5f,.5f,1.4f}}};
    auto adjusted=zTsParametricStairs::optimiseGraph(levelGraph,os,.25,.35,.15,.18);zFnGraph(adjusted.graph).getVertexPositions(levels);
    check(levels[0].z==0&&std::abs(levels[1].z-.9)<1e-6&&adjusted.adjustedLandings==1,"Bottom anchored while infeasible landing moves to nearest range-feasible level");
    check(adjusted.riserCounts==zIntArray{5}&&std::abs(adjusted.treadDepths[0]-.25)<1e-6&&std::abs(adjusted.riserHeights[0]-.18)<1e-6,"Optimiser selects feasible flight tread and riser within domains");
    check(std::abs(adjusted.landingPolygons[0][0].z-.9)<1e-6,"Imported landing outline follows optimised level");
    os.landingPolygons=adjusted.landingPolygons;os.landingVertices=adjusted.landingVertices;os.flightRiserCounts=adjusted.riserCounts;os.flightTreadDepths=adjusted.treadDepths;t.update(adjusted.graph,os);check(t.flightRiserCounts()==adjusted.riserCounts,"Selected dimensions drive generated geometry");
    rejects([&]{zTsParametricStairs::optimiseGraph(levelGraph,os,.4,.2,.15,.18);});
    zObjectGraph valley;zPointArray vp={{0,0,1},{4,0,0},{8,0,1}};zIntArray ve={0,1,1,2};zFnGraph(valley).create(vp,ve);os.landingPolygons.clear();os.landingVertices.clear();os.flightRiserCounts.clear();os.flightTreadDepths.clear();
    auto valleyPlan=zTsParametricStairs::optimiseGraph(valley,os,.25,.35,.15,.18);check(valleyPlan.bottomNodes==zIntArray{1},"Lowest degree-two node becomes anchored landing");zFnGraph(valleyPlan.graph).getVertexPositions(levels);check(levels[1].z==0,"Valley bottom elevation remains fixed");
    }
    std::cout<<"Parametric Stairs native checks passed\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
