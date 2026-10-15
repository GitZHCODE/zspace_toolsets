#include <zspace/toolsets.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace zSpace;
void check(bool b,const char* text){if(!b)throw std::runtime_error(text);}
template<class F> void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"Expected rejection");}
zObjectMesh arch(float z=1){zObjectMesh m;zPointArray p={{-1,0,0},{1,0,0},{0,0,z}};zIntArray c={3},f={0,1,2};zFnMesh(m).create(p,c,f);return m;}
int main(){try{
    zObjectMesh unequal;zPointArray up={{-1,0,0},{1,0,0},{0.6f,0,1}};
    zIntArray uc={3},uf={0,1,2};zFnMesh(unequal).create(up,uc,uf);
    double freeDrift=0;
    for(int mode=0;mode<3;++mode){
        zTsVault plan;plan.setInputMesh(unequal);plan.setSupports({0,1});plan.setVertexMasses({0,0,2});
        zVaultBestFitSettings ps;ps.planConstraint=static_cast<zVaultPlanConstraint>(mode);ps.xyPenaltyWeight=20;
        ps.minMaxEdgeRatio=.9;ps.edgeRegularizationStrength=1;plan.setBestFitSettings(ps);
        plan.computeBestFitForceDensities();auto r=plan.stepBestFit(20);
        if(mode==0){freeDrift=r.rmsXYDisplacement;check(freeDrift>0,"Free plan moves horizontally");}
        if(mode==1)check(r.rmsXYDisplacement>0&&r.rmsXYDisplacement<freeDrift,"Soft XY reduces drift without locking plan");
        if(mode==2)check(r.rmsXYDisplacement==0,"Plan mode hard lock");
    }
    zTsVault reference;reference.setInputMesh(unequal);reference.setSupports({0,1});reference.setVertexMasses({0,0,2});
    auto target=up;target[2].x=.2f;reference.setReferencePositions(target);
    zVaultBestFitSettings ps;ps.planConstraint=zVaultPlanConstraint::PenalizeXY;ps.xyPenaltyWeight=2;reference.setBestFitSettings(ps);
    auto refStart=reference.computeBestFitForceDensities();check(refStart.rmsXYDisplacement>.39&&!refStart.converged,"Separate reference evaluated at initial fit");
    auto refEnd=reference.stepBestFit(10);check(refEnd.rmsXYDisplacement<refStart.rmsXYDisplacement,"Attraction towards persistent reference");
    reference.previousBestFit();reference.resetBestFit();check(reference.report().rmsXYDisplacement==refStart.rmsXYDisplacement,"Reset does not rebase reference");
    ps.xyPenaltyWeight=3;reference.setBestFitSettings(ps);check(reference.computeBestFitForceDensities().rmsXYDisplacement==refStart.rmsXYDisplacement,"Weight edits retain reference");
    rejects([&]{reference.setReferencePositions({zPoint()});});ps.xyPenaltyWeight=-1;rejects([&]{reference.setBestFitSettings(ps);});
    for(bool fix:{false,true}){
        zTsVault regular;regular.setInputMesh(unequal);regular.setSupports({0,1});regular.setVertexMasses({0,0,2});
        zVaultBestFitSettings rs;rs.fixXY=fix;rs.minMaxEdgeRatio=.9;rs.edgeRegularizationStrength=1;
        regular.setBestFitSettings(rs);auto initial=regular.computeBestFitForceDensities();
        check(!initial.converged&&initial.residualNorm<1e-6,"Regularization runs on a physically balanced but unequal mesh");
        auto after=regular.stepBestFit(20);
        check(after.iteration>0&&after.minMaxEdgeRatio>initial.minMaxEdgeRatio,"Regularization improves edge ratio");
        check(after.edgeRegularizationError<initial.edgeRegularizationError,"Regularization decreases length penalty");
        zPointArray p;zFnMesh(regular.resultMesh()).getVertexPositions(p);
        check(p[0].x==-1&&p[1].x==1&&p[0].z==0&&p[1].z==0,"Regularizer preserves supports");
        if(fix)check(p[2].x==up[2].x&&p[2].y==up[2].y&&p[2].z!=up[2].z,"Regularizer respects Fix XY");
        check(regular.previousBestFit(),"Regularization history");regular.resetBestFit();
        check(std::abs(regular.report().edgeRegularizationError-initial.edgeRegularizationError)<1e-12,"Reset restores length error");
        rs.edgeRegularizationStrength=0;regular.setBestFitSettings(rs);check(regular.computeBestFitForceDensities().converged,"Zero strength preserves equilibrium-only behaviour");
        rs.minMaxEdgeRatio=1.1;rejects([&]{regular.setBestFitSettings(rs);});rs.minMaxEdgeRatio=.5;rs.edgeRegularizationStrength=-1;rejects([&]{regular.setBestFitSettings(rs);});
    }
    // Asymmetric free vertex has horizontal residual after q is bounded.
    // The unconstrained control must move in plan, while Fix XY must not.
    zObjectMesh asymmetric;zPointArray ap={{-1,-1,0},{1,1,0},{0.4f,0.2f,1}};
    zIntArray ac={3},af={0,1,2};zFnMesh(asymmetric).create(ap,ac,af);
    for(bool fix:{false,true}){
        zTsVault constrained;constrained.setInputMesh(asymmetric);constrained.setSupports({0,1});constrained.setVertexMasses({0,0,2});
        zVaultBestFitSettings s;s.fixXY=fix;s.qMin=0.1;s.qMax=0.2;s.relativeTolerance=1e-6;
        constrained.setBestFitSettings(s);auto first=constrained.computeBestFitForceDensities();
        auto next=constrained.stepBestFit(10);zPointArray points;zFnMesh(constrained.resultMesh()).getVertexPositions(points);
        check(next.iteration>0&&next.residualNorm<first.residualNorm,"Constrained solve makes progress");
        check(points[2].z!=ap[2].z,"Free Z updates");
        check(points[0].x==ap[0].x&&points[0].y==ap[0].y&&points[0].z==ap[0].z&&points[1].z==ap[1].z,"Supports stay fixed");
        if(fix){
            check(points[2].x==ap[2].x&&points[2].y==ap[2].y,"Fix XY preserves exact input coordinates");
            check(constrained.previousBestFit(),"Constrained rewind");zFnMesh(constrained.resultMesh()).getVertexPositions(points);
            check(points[2].x==ap[2].x&&points[2].y==ap[2].y,"Rewind preserves XY");
            constrained.resetBestFit();zFnMesh(constrained.resultMesh()).getVertexPositions(points);
            check(points[2].z==ap[2].z&&points[2].x==ap[2].x&&points[2].y==ap[2].y,"Constrained reset");
        }else check(points[2].x!=ap[2].x&&points[2].y!=ap[2].y,"Default mode can move XY");
    }
    auto mesh=arch();zTsVault v;v.setInputMesh(mesh);v.setSupports({0,1});v.setVertexMasses({0,0,2});
    auto r=v.computeBestFitForceDensities();check(r.converged,"Analytic compression arch equilibrium");
    for(int e:v.activeEdgeIds())check(std::abs(v.forceDensities()[e]+1)<1e-8,"q=-1 for mass 2, gravity -1");
    check(std::abs(v.supportReactions()[0].z-1)<1e-6,"Support reaction");
    check(v.activeEdgeIds().size()==2,"Support-support edge excluded");
    auto settings=zVaultBestFitSettings();settings.qMin=0.1;settings.qMax=0.2;settings.relativeTolerance=1e-3;v.setBestFitSettings(settings);
    r=v.computeBestFitForceDensities();double start=r.residualNorm;
    auto initial=v.forceDensities();r=v.stepBestFit(10);check(r.residualNorm<start,"Compression iteration reduces residual");
    zPointArray p;zFnMesh(v.resultMesh()).getVertexPositions(p);check(p[0].x==-1 && p[0].z==0 && p[1].x==1 && p[1].z==0,"Supports preserved");
    for(int e:v.activeEdgeIds())check(v.forceDensities()[e]>=-0.200000001&&v.forceDensities()[e]<=-0.099999999,"Bounds respected");
    int iteration=r.iteration;check(v.previousBestFit()&&v.report().iteration==iteration-1,"Rewind accepted state");
    v.resetBestFit();check(v.report().iteration==0&&std::abs(v.report().residualNorm-start)<1e-10&&v.forceDensities()==initial,"Reset geometry and densities");
    rejects([&]{v.setVertexMass(-1);});rejects([&]{v.setSupports({99});});rejects([&]{v.setVertexMasses({1});});
    settings.qMin=2;settings.qMax=1;rejects([&]{v.setBestFitSettings(settings);});
    auto cable=arch(-1);v.setInputMesh(cable);v.setSupports({0,1});v.setVertexMasses({0,0,2});
    zIntArray tension;auto e=v.edgeConnects();for(int i=0;i<int(e.size()/2);++i)if(e[2*i]==2||e[2*i+1]==2)tension.push_back(i);
    v.setTensionEdges(tension);r=v.computeBestFitForceDensities();check(r.converged,"Tension cable equilibrium");
    for(int edge:tension)check(std::abs(v.forceDensities()[edge]-1)<1e-8,"Tension sign");
    v.setVertexMass(0);rejects([&]{v.computeBestFitForceDensities();});
    // A planar, loaded mechanism has rank deficiency: fit remains finite, no false equilibrium.
    auto flat=arch(0);v.setInputMesh(flat);v.setSupports({0,1});r=v.computeBestFitForceDensities();
    check(r.rank<r.activeEdges&&!r.converged,"Rank deficiency diagnostic");
    for(double q:v.forceDensities())check(std::isfinite(q),"Finite rank-deficient fit");
    // Mixed-sign star: lower supports push up, upper supports pull up.
    zObjectMesh mixed;zPointArray mp={{-1,0,-1},{0,-1,1},{1,0,-1},{0,1,1},{0,0,0}};
    zIntArray mc={3,3,3,3},mf={0,1,4,1,2,4,2,3,4,3,0,4};zFnMesh(mixed).create(mp,mc,mf);
    v.setInputMesh(mixed);v.setSupports({0,1,2,3});v.setVertexMasses({0,0,0,0,4});
    tension.clear();e=v.edgeConnects();for(int i=0;i<int(e.size()/2);++i)
        if((e[2*i]==4&&(e[2*i+1]==1||e[2*i+1]==3))||(e[2*i+1]==4&&(e[2*i]==1||e[2*i]==3)))tension.push_back(i);
    v.setTensionEdges(tension);r=v.computeBestFitForceDensities();check(r.converged,"Mixed tension/compression equilibrium");
    for(int edge:v.activeEdgeIds())check(std::abs(std::abs(v.forceDensities()[edge])-1)<1e-8,"Mixed forces known magnitude");
    auto savedQ=v.forceDensities();
    auto diagram=v.createForceFromForm();check(diagram.exact,"Mixed-sign reciprocal closes");
    check(savedQ==v.forceDensities(),"Diagram reconstruction does not change q");
    zObjectMesh form=v.formMesh(),force=v.forceMesh();
    for(zItMeshHalfEdge h(form);!h.end();h++){
        int f=v.formHalfEdgeToForceHalfEdge()[h.getId()];if(f<0)continue;
        check(v.forceHalfEdgeToFormHalfEdge()[f]==h.getId(),"Reciprocal halfedge mapping round trip");
        zItMeshHalfEdge dual(force,f);auto d=h.getVector(),fd=dual.getVector();
        int edge=-1;zFnMesh(v.resultMesh()).edgeExists(h.getStartVertex().getId(),h.getVertex().getId(),edge);
        check(std::abs((fd*d)/(d*d)-savedQ[edge])<1e-6,"Signed q recovered from reciprocal edge vectors");
    }
    zDoubleArray lower(e.size()/2,-2),upper(e.size()/2,2);lower[tension[0]]=0.5;upper[tension[0]]=0.75;
    v.setForceDensityBounds(lower,upper);v.computeBestFitForceDensities();check(v.forceDensities()[tension[0]]==0.75,"Per-edge signed bound applied");
    diagram=v.updateForceFromDensities();check(!diagram.exact&&diagram.maxEdgeMismatch>0&&diagram.closureError>0,"Unbalanced fit is explicitly approximate");
    upper[tension[0]]=-0.1;lower[tension[0]]=-1;v.setForceDensityBounds(lower,upper);rejects([&]{v.computeBestFitForceDensities();});
    // Free boundary vertices share exterior regions only between supports.
    zObjectMesh grid;zPointArray gp;for(int y=0;y<3;++y)for(int x=0;x<3;++x)gp.emplace_back(float(x),float(y),float(x==1&&y==1));
    zIntArray gc={4,4,4,4},gf={0,1,4,3,1,2,5,4,3,4,7,6,4,5,8,7};zFnMesh(grid).create(gp,gc,gf);
    zTsVault boundary;boundary.setInputMesh(grid);boundary.setSupports({0,2,6,8});
    rejects([&]{boundary.createForceFromForm();});boundary.computeBestFitForceDensities();
    auto br=boundary.createForceFromForm();check(std::isfinite(br.maxEdgeMismatch),"Supported boundary segment reconstruction");
    zObjectMesh bf=boundary.forceMesh();check(zFnMesh(bf).numVertices()==8,"Four source faces plus four exterior regions");
    auto colors=zColorArray();zFnMesh(bf).getEdgeColors(colors);check(!colors.empty()&&std::abs(colors[0].b-240.0f/255)<1e-6,"Compression diagram edge colors");
    boundary.setGravity(zVector(1,0,-1));boundary.computeBestFitForceDensities();rejects([&]{boundary.createForceFromForm();});
    // An unsupported inner opening is a valid exterior dual region.
    zObjectMesh ring;zPointArray rp={{-2,-2,0},{2,-2,0},{2,2,0},{-2,2,0},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    zIntArray rc={4,4,4,4},rf={0,1,5,4,1,2,6,5,2,3,7,6,3,0,4,7};zFnMesh(ring).create(rp,rc,rf);
    zTsVault opening;opening.setInputMesh(ring);opening.setSupports({0,1,2,3});
    check(opening.computeBestFitForceDensities().converged,"Annular vault equilibrium");
    auto oq=opening.forceDensities();check(opening.createForceFromForm().exact,"Unsupported hole reciprocal closes");
    check(opening.forceDensities()==oq,"Opening reconstruction preserves fitted q");
    zObjectMesh of=opening.forceMesh();check(zFnMesh(of).numVertices()==5,"Hole produces one dual vertex, not artificial supports");
    int mapped=0;for(int id:opening.formHalfEdgeToForceHalfEdge())if(id>=0)++mapped;
    check(mapped==2*opening.activeEdgeIds().size(),"Every active hole/radial edge has a reciprocal");
    // Signed axial sums cancel, so the original geometry direction is exactly zero.
    // Fixed q bounds also eliminate its q direction, but equilibrium exists at z=-1.
    for(double weight:{0.5,1.0}){
        zObjectMesh stalled;zPointArray sp={{-2,0,0},{0,-1,0},{2,0,0},{0,1,0},{0,0,0}};
        zFnMesh(stalled).create(sp,mc,mf);zTsVault rescue;rescue.setInputMesh(stalled);rescue.setSupports({0,1,2,3});
        zDoubleArray bounds(rescue.edgeConnects().size()/2,-0.5);zIntArray tensile;
        auto se=rescue.edgeConnects();for(int j=0;j<int(bounds.size());++j)
            if((se[2*j]==4&&(se[2*j+1]==1||se[2*j+1]==3))||(se[2*j+1]==4&&(se[2*j]==1||se[2*j]==3))){bounds[j]=1;tensile.push_back(j);}
        rescue.setTensionEdges(tensile);rescue.setForceDensityBounds(bounds,bounds);
        auto config=zVaultBestFitSettings();config.geometryWeight=weight;config.fixXY=true;rescue.setBestFitSettings(config);
        auto initialReport=rescue.computeBestFitForceDensities();check(initialReport.residualNorm>0.9,"Cancellation fixture starts unbalanced");
        auto solved=rescue.stepBestFit(1);check(solved.converged&&solved.iteration==1,"Fallback escapes cancelled mixed-sign direction");
        zPointArray finalPoints;zFnMesh(rescue.resultMesh()).getVertexPositions(finalPoints);
        for(int i=0;i<4;++i)check((finalPoints[i]-sp[i]).length()<1e-6,"Fallback fixes supports");
        check(finalPoints[4].x==0&&finalPoints[4].y==0,"Fallback fixes XY");
        check(rescue.forceDensities()[tensile[0]]==1,"Fallback respects fixed q bounds");
        check(rescue.previousBestFit()&&rescue.report().iteration==0,"Fallback can rewind");
        rescue.stepBestFit(1);rescue.resetBestFit();check(rescue.report().iteration==0&&rescue.report().residualNorm>0.9,"Fallback can reset");
    }
    std::cout<<"Vault analytic fit, diagrams, signed reciprocal round trip, boundaries, history and validation passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
