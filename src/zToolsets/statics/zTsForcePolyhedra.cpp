#include "zTsForcePolyhedraInternal.h"
#include <algorithm>
#include <cmath>
#include <set>

namespace zSpace {
namespace forcePolyhedra {
void alignVirtualCaps(Cell& cell, bool closeTriangle);
}
using namespace forcePolyhedra;
namespace {
double degrees(const V& a,const V& b){return std::atan2(a.cross(b).norm(),a.dot(b))*180.0/3.14159265358979323846;}
V blend(const V& a,const V& b,double weight){V n=(1-weight)*a+weight*b;return n.norm()>1e-12?n.normalized():a;}

}
struct zTsForcePolyhedra::Impl {
    zPointArray original,positions;
    std::vector<V> precise;
    zIntArray edges,supports,tension;
    zDoubleArray masses,forces,magnitudes;
    zVector gravity=zVector(0,0,-1);
    zVectorArray loads,residuals,reactions;
    zObjectGraph graph;
    std::unique_ptr<GraphDynamicsState> formDynamics;
    zForcePolyhedraSettings settings;
    zForcePolyhedraReport report;
    std::vector<Cell> cells;
    zIntArray memberFaces; // Two indices per edge, in its two node cells.
    double maximumLength=1;
    void invalidate(){report={};cells.clear();memberFaces.clear();forces.clear();residuals.clear();reactions.clear();formDynamics.reset();}
    void ready()const{require(!original.empty(),"Set an input graph first.");}
    std::vector<bool> fixed()const{std::vector<bool> f(positions.size(),false);for(int i:supports)f[i]=true;return f;}
    Eigen::MatrixXd external()const{
        Eigen::MatrixXd p(positions.size(),3);
        for(int i=0;i<int(positions.size());++i)p.row(i)=(vector(loads[i])+masses[i]*vector(gravity)).transpose();
        require(p.allFinite(),"Nonfinite total load.");return p;
    }
    void updateGraph(){for(int i=0;i<int(precise.size());++i)positions[i]=publicVector(precise[i]);auto p=positions;auto e=edges;zFnGraph(graph).create(p,e);}
    int checkedCell(int id)const{require(id>=0&&id<int(cells.size()),"Force cell index out of range; create force polyhedra first.");return id;}
    double memberArea(int e,int side)const{
        int node=edges[2*e+side],index=memberFaces[2*e+side];const auto& f=cells[node].faces[index];
        double total=0,share=0;for(auto source:f.sources){total+=source.magnitude;if(source.type==zForcePolyhedraFaceType::Member&&source.id==e)share=source.magnitude;}
        return total>0?f.area*share/total:0;
    }
    void summarize(){
        auto p=external();auto pinned=fixed();report.maxResidual=report.residualNorm=0;
        report.maxAngleDeviation=report.maxContactAreaMismatch=0;
        report.maxCorrespondingAreaDifference=0;report.worstAreaDifferenceEdge=-1;
        report.cellsConverged=true;report.maxAreaError=report.maxNormalError=report.maxPlanarityError=0;report.iterations=0;
        report.maxCellVertexStep=0;report.worstAreaCell=report.worstAngleEdge=report.worstMotionCell=-1;
        for(int cellId=0;cellId<int(cells.size());++cellId){auto& c=cells[cellId];
            refresh(c,settings);report.cellsConverged=report.cellsConverged&&c.report.converged;
            if(c.report.maxAreaError>report.maxAreaError){report.maxAreaError=c.report.maxAreaError;report.worstAreaCell=cellId;}
            if(c.report.maxVertexStep>report.maxCellVertexStep){report.maxCellVertexStep=c.report.maxVertexStep;report.worstMotionCell=cellId;}
            report.maxNormalError=std::max(report.maxNormalError,c.report.maxNormalError);
            report.maxPlanarityError=std::max(report.maxPlanarityError,c.report.maxPlanarityError);
            report.iterations=std::max(report.iterations,c.report.iterations);
        }
        for(int e=0;e<int(forces.size());++e){
            int u=edges[2*e],v=edges[2*e+1];V edge=precise[v]-precise[u];
            require(edge.norm()>1e-10,"Collapsed graph member.");
            V force=forces[e]*edge.normalized();p.row(u)+=force.transpose();p.row(v)-=force.transpose();
            V n0=faceNormal(cells[u],memberFaces[2*e]),n1=faceNormal(cells[v],memberFaces[2*e+1]);
            V direction=edge.normalized()*(forces[e]>0?1.0:-1.0);
            double error=std::max({degrees(direction,n0),degrees(-direction,n1),degrees(n0,-n1)});
            if(error>report.maxAngleDeviation){report.maxAngleDeviation=error;report.worstAngleEdge=e;}
            double a=memberArea(e,0),b=memberArea(e,1),scale=std::max(std::abs(forces[e])*settings.forceScale,1e-15);
            report.maxContactAreaMismatch=std::max(report.maxContactAreaMismatch,std::abs(a-b)/scale);
            if(std::abs(a-b)>report.maxCorrespondingAreaDifference){report.maxCorrespondingAreaDifference=std::abs(a-b);report.worstAreaDifferenceEdge=e;}
        }
        require(p.allFinite(),"Nonfinite equilibrium result.");residuals.assign(positions.size(),zVector());reactions.assign(positions.size(),zVector());
        for(int i=0;i<int(positions.size());++i){V r=p.row(i).transpose();if(pinned[i])reactions[i]=publicVector(-r);else{
            residuals[i]=publicVector(r);report.maxResidual=std::max(report.maxResidual,r.norm());report.residualNorm+=r.squaredNorm();}}
        report.residualNorm=std::sqrt(report.residualNorm);
        report.areaPursuitResolved=report.maxCorrespondingAreaDifference<=settings.correspondingAreaTolerance;
        report.equilibrium=report.cellsConverged&&report.maxResidual<=settings.equilibriumTolerance&&
            report.maxAngleDeviation<=settings.angleTolerance&&report.areaPursuitResolved;
        report.status=report.equilibrium?"equilibrium":"reciprocation_incomplete";
    }
    // Legacy zTsPolytopal: targets are oriented face normals; each free vertex
    // approaches the average of neighbour + target * clamped current length.
    // Residuals are integrated by zFnGraphDynamics, retaining exact pins.
    void moveForm(const std::vector<V>& directions){
        report.maxFormVertexStep=0;
        if(settings.formWeight==0)return;
        if(!formDynamics){
            formDynamics=std::make_unique<GraphDynamicsState>();auto& state=*formDynamics;state.scale=maximumLength;
            zPointArray points;for(const V& p:precise)points.push_back(publicVector((p-precise[0])/state.scale));auto connects=edges;
            zFnGraph(state.object).create(points,connects);state.fn.makeDynamic(false);auto pins=supports;state.fn.setFixed(pins);
        }
        auto pinned=fixed();std::vector<V> sums(positions.size(),V::Zero());std::vector<int> count(positions.size(),0);
        for(int e=0;e<int(forces.size());++e){
            int u=edges[2*e],v=edges[2*e+1];V a=precise[u],b=precise[v];
            double length=std::clamp((b-a).norm(),settings.minMaxEdgeRatio*maximumLength,maximumLength);
            sums[u]+=b-directions[e]*length;sums[v]+=a+directions[e]*length;++count[u];++count[v];
        }
        std::vector<V> corrections(precise.size(),V::Zero());
        for(int i=0;i<int(precise.size());++i)if(!pinned[i]&&count[i])corrections[i]=sums[i]/count[i]-precise[i];
        auto delta=integrateCorrections(*formDynamics,corrections,settings,settings.maxVertexStepRatio*maximumLength);auto next=precise;
        for(int i=0;i<int(next.size());++i)if(!pinned[i]){next[i]+=delta[i];report.maxFormVertexStep=std::max(report.maxFormVertexStep,delta[i].norm());}
        for(int e=0;e<int(forces.size());++e)if((next[edges[2*e]]-next[edges[2*e+1]]).norm()<=1e-10)return;
        precise=std::move(next);
    }
};
zTsForcePolyhedra::zTsForcePolyhedra():impl(std::make_unique<Impl>()){}
zTsForcePolyhedra::~zTsForcePolyhedra()=default;
void zTsForcePolyhedra::clear(){impl=std::make_unique<Impl>();}
void zTsForcePolyhedra::setInputGraph(zObjectGraph& graph){
    auto next=std::make_unique<Impl>();zFnGraph fn(graph);fn.getVertexPositions(next->original);fn.getEdgeData(next->edges);
    require(!next->original.empty()&&next->original.size()<=1024&&next->edges.size()<=4096,"Expected 1..1024 vertices and at most 2048 edges.");
    for(auto p:next->original)require(vector(p).allFinite(),"Nonfinite graph coordinate.");
    require(next->edges.size()%2==0,"Invalid graph edge array.");std::set<std::pair<int,int>> seen;next->maximumLength=0;
    for(size_t e=0;e<next->edges.size();e+=2){int u=next->edges[e],v=next->edges[e+1];
        require(u>=0&&v>=0&&u<int(next->original.size())&&v<int(next->original.size())&&u!=v,"Invalid graph member.");
        require(seen.emplace(std::min(u,v),std::max(u,v)).second,"Duplicate graph member.");
        double length=(vector(next->original[u])-vector(next->original[v])).norm();require(length>1e-10,"Collapsed graph member.");
        next->maximumLength=std::max(next->maximumLength,length);
    }
    next->positions=next->original;for(auto p:next->positions)next->precise.push_back(vector(p));next->masses.assign(next->positions.size(),0);next->loads.assign(next->positions.size(),zVector());next->updateGraph();impl=std::move(next);
}
void zTsForcePolyhedra::setSupports(const zIntArray& ids){
    impl->ready();std::set<int> seen;for(int id:ids)require(id>=0&&id<int(impl->positions.size())&&seen.insert(id).second,"Invalid or duplicate support ID.");
    impl->supports=ids;impl->invalidate();
}
void zTsForcePolyhedra::setFormPositions(const zPointArray& points){
    auto& s=*impl;s.ready();require(points.size()==s.positions.size(),"Expected one form position per vertex.");auto pinned=s.fixed();
    for(int i=0;i<int(points.size());++i){require(vector(points[i]).allFinite(),"Nonfinite form coordinate.");if(pinned[i])require((vector(points[i])-vector(s.original[i])).norm()==0,"Pinned form vertices cannot move.");}
    for(size_t e=0;e<s.edges.size();e+=2)require((vector(points[s.edges[e]])-vector(points[s.edges[e+1]])).norm()>1e-10,"Collapsed form edge.");
    s.positions=points;s.precise.clear();for(auto p:points)s.precise.push_back(vector(p));s.formDynamics.reset();s.updateGraph();if(!s.cells.empty())s.summarize();
}
void zTsForcePolyhedra::setTensionEdges(const zIntArray& ids){
    impl->ready();std::set<int> seen;for(int id:ids)require(id>=0&&id<int(impl->edges.size()/2)&&seen.insert(id).second,"Invalid or duplicate tension edge ID.");
    impl->tension=ids;impl->invalidate();
}
void zTsForcePolyhedra::setEdgeForceMagnitudes(const zDoubleArray& values){
    impl->ready();require(values.empty()||values.size()==impl->edges.size()/2,"Expected one target force magnitude per edge.");
    for(double f:values)require(std::isfinite(f)&&f>0&&f<=1e20,"Force magnitudes must be positive and finite.");impl->magnitudes=values;impl->invalidate();
}
void zTsForcePolyhedra::setVertexMass(double mass){impl->ready();require(std::isfinite(mass)&&mass>=0,"Mass must be finite and nonnegative.");setVertexMasses(zDoubleArray(impl->positions.size(),mass));}
void zTsForcePolyhedra::setVertexMasses(const zDoubleArray& masses){
    impl->ready();require(masses.size()==impl->positions.size(),"Expected one mass per vertex.");for(double x:masses)require(std::isfinite(x)&&x>=0,"Mass must be finite and nonnegative.");impl->masses=masses;impl->invalidate();
}
void zTsForcePolyhedra::setGravity(const zVector& g){require(vector(g).allFinite(),"Gravity must be finite.");impl->gravity=g;impl->invalidate();}
void zTsForcePolyhedra::setVertexLoads(const zVectorArray& loads){
    impl->ready();require(loads.size()==impl->positions.size(),"Expected one load per vertex.");for(const auto& p:loads)require(vector(p).allFinite(),"Loads must be finite.");impl->loads=loads;impl->invalidate();
}
void zTsForcePolyhedra::setSettings(const zForcePolyhedraSettings& s){
    require(std::isfinite(s.equilibriumTolerance)&&s.equilibriumTolerance>0&&std::isfinite(s.areaTolerance)&&s.areaTolerance>0&&
        std::isfinite(s.correspondingAreaTolerance)&&s.correspondingAreaTolerance>=0&&
        std::isfinite(s.normalTolerance)&&s.normalTolerance>0&&s.normalTolerance<1&&std::isfinite(s.forceScale)&&s.forceScale>0&&
        s.maxIterations>0&&s.maxIterations<=100000&&std::isfinite(s.formWeight)&&s.formWeight>=0&&s.formWeight<=1&&
        std::isfinite(s.angleTolerance)&&s.angleTolerance>0&&s.angleTolerance<90&&std::isfinite(s.relaxation)&&s.relaxation>0&&s.relaxation<=1&&
        std::isfinite(s.minMaxEdgeRatio)&&s.minMaxEdgeRatio>0&&s.minMaxEdgeRatio<=1&&s.areaIterationsPerStep>0&&s.areaIterationsPerStep<=1000&&
        std::isfinite(s.timeStep)&&s.timeStep>0&&s.timeStep<=1&&std::isfinite(s.maxVertexStepRatio)&&s.maxVertexStepRatio>0&&s.maxVertexStepRatio<=.25&&
        std::isfinite(s.normalStepDegrees)&&s.normalStepDegrees>0&&s.normalStepDegrees<=30,
        "Invalid force polyhedra settings.");
    // Movement-weight edits can continue an existing form/force session.
    bool reconstruct=s.forceScale!=impl->settings.forceScale;
    impl->settings=s;if(reconstruct)impl->invalidate();else if(!impl->cells.empty())impl->summarize();
}
zForcePolyhedraReport zTsForcePolyhedra::createForcePolyhedra(){
    auto& s=*impl;s.ready();auto p=s.external(),net=p;const auto pinned=s.fixed();
    std::vector<Cell> next(s.positions.size());zDoubleArray forces(s.edges.size()/2);
    double seed=1,totalLoad=0;for(int i=0;i<p.rows();++i){seed=std::max(seed,p.row(i).norm());totalLoad+=p.row(i).norm();}
    // A long loaded arch carries accumulated nodal loads towards its supports.
    // Seeding every member from only one node's weight produces large initial
    // normal rotations and an unnecessarily poor reciprocal starting diagram.
    seed=std::max(seed,totalLoad/std::max(size_t(1),s.supports.size()));
    std::set<int> tension(s.tension.begin(),s.tension.end());
    for(int e=0;e<int(forces.size());++e)forces[e]=(tension.count(e)?1.0:-1.0)*(s.magnitudes.empty()?seed:s.magnitudes[e]);
    auto add=[&](int node,const V& force,zForcePolyhedraFaceType type,int id){
        double magnitude=force.norm();if(magnitude<1e-12)return;auto& c=next[node];V normal=force/magnitude;int index=-1;
        for(int j=0;j<int(c.normals.size());++j)if((c.normals[j]-normal).norm()<1e-9){index=j;break;}
        if(index<0){index=int(c.faces.size());c.normals.push_back(normal);c.faces.push_back({});c.faces.back().normal=publicVector(normal);c.faces.back().type=type;}
        c.faces[index].targetArea+=magnitude*s.settings.forceScale;c.faces[index].sources.push_back({type,id,magnitude});
    };
    for(int e=0;e<int(forces.size());++e){int u=s.edges[2*e],v=s.edges[2*e+1];V direction=s.precise[v]-s.precise[u];V force=direction.normalized()*forces[e];
        add(u,force,zForcePolyhedraFaceType::Member,e);add(v,-force,zForcePolyhedraFaceType::Member,e);net.row(u)+=force.transpose();net.row(v)-=force.transpose();}
    for(int i=0;i<int(next.size());++i){add(i,p.row(i).transpose(),zForcePolyhedraFaceType::Load,i);if(pinned[i])add(i,-net.row(i).transpose(),zForcePolyhedraFaceType::Reaction,i);}
    for(int i=0;i<int(next.size());++i){auto& c=next[i];try{buildCell(c,s.settings,true);pursue(c,s.settings,s.settings.maxIterations,c.report.closureError>s.settings.equilibriumTolerance);}catch(const std::exception& e){throw std::runtime_error("Force cell at graph vertex "+std::to_string(i)+": "+e.what());}}
    zIntArray mappings(2*forces.size(),-1);
    for(int i=0;i<int(next.size());++i)for(int j=0;j<int(next[i].faces.size());++j)for(auto source:next[i].faces[j].sources)if(source.type==zForcePolyhedraFaceType::Member){int e=source.id;int side=s.edges[2*e]==i?0:1;mappings[2*e+side]=j;}
    for(int id:mappings)require(id>=0,"Missing member-to-cell face mapping.");
    s.cells=std::move(next);s.memberFaces=std::move(mappings);s.forces=std::move(forces);s.report={};s.summarize();return s.report;
}
void zTsForcePolyhedra::setFormWeight(double weight){auto settings=impl->settings;settings.formWeight=weight;setSettings(settings);}
zForcePolyhedraReport zTsForcePolyhedra::stepAreaPursuit(int iterations){
    require(iterations>0&&iterations<=100000,"Invalid area-pursuit iteration count.");require(!impl->cells.empty(),"Create force polyhedra first.");
    for(auto& c:impl->cells)pursue(c,impl->settings,iterations,c.report.closureError>impl->settings.equilibriumTolerance);impl->summarize();return impl->report;
}
zForcePolyhedraReport zTsForcePolyhedra::stepEquilibrium(int iterations){
    require(iterations>0&&iterations<=100000,"Invalid equilibrium iteration count.");auto& s=*impl;require(!s.cells.empty(),"Create force polyhedra before reciprocation.");
    const int m=int(s.forces.size());
    for(int iteration=0;iteration<iterations&&!s.report.equilibrium;++iteration){
        std::vector<V> directions(m);std::vector<std::vector<V>> proposals(s.cells.size());
        std::vector<std::vector<int>> counts(s.cells.size());
        for(int i=0;i<int(s.cells.size());++i){proposals[i].assign(s.cells[i].faces.size(),V::Zero());counts[i].assign(s.cells[i].faces.size(),0);}
        for(int e=0;e<m;++e){
            int u=s.edges[2*e],v=s.edges[2*e+1],f0=s.memberFaces[2*e],f1=s.memberFaces[2*e+1];
            V n0=faceNormal(s.cells[u],f0),n1=-faceNormal(s.cells[v],f1);
            V common=blend(n0,n1,.5);double sign=s.forces[e]>0?1.0:-1.0;
            V graphNormal=(s.precise[v]-s.precise[u]).normalized()*sign;
            // formWeight controls how much each diagram moves, independently of
            // the legacy interpolation convention: endpoints freeze one diagram.
            directions[e]=blend(graphNormal,common,s.settings.formWeight)*sign;
            V forceTarget=blend(common,graphNormal,1-s.settings.formWeight);
            proposals[u][f0]+=forceTarget;proposals[v][f1]-=forceTarget;++counts[u][f0];++counts[v][f1];
        }
        if(s.settings.formWeight<1){
            if(s.magnitudes.empty()){
                // Unspecified member forces are free to redistribute. Use the
                // mean of the two actual contact contributions as their common
                // next area target; prescribed nodal load areas remain fixed.
                zDoubleArray next(m);for(int e=0;e<m;++e){double target=(s.memberArea(e,0)+s.memberArea(e,1))/(2*s.settings.forceScale);next[e]=std::max(1e-10,(1-s.settings.relaxation)*std::abs(s.forces[e])+s.settings.relaxation*target);}
                for(int e=0;e<m;++e)s.forces[e]=(s.forces[e]>0?1:-1)*next[e];
                for(auto& c:s.cells)for(auto& f:c.faces)if(!f.sources.empty()){
                    f.targetArea=0;for(auto& source:f.sources){if(source.type==zForcePolyhedraFaceType::Member)source.magnitude=next[source.id];f.targetArea+=source.magnitude*s.settings.forceScale;}
                }
            }
            for(int i=0;i<int(s.cells.size());++i){auto& c=s.cells[i];
                for(int f=0;f<int(c.faces.size());++f){
                    bool load=false,reaction=false;for(auto source:c.faces[f].sources){load|=source.type==zForcePolyhedraFaceType::Load;reaction|=source.type==zForcePolyhedraFaceType::Reaction;}
                    if(load)continue; // Preserve applied/self-weight load directions.
                    if(counts[i][f]&&proposals[i][f].norm()>1e-12)c.normals[f]=proposals[i][f].normalized();
                    else if(reaction)c.normals[f]=faceNormal(c,f);
                    c.faces[f].normal=publicVector(c.normals[f]);
                }
                alignVirtualCaps(c,false);
                // At a pinned node the reaction is the closing vector, not an
                // independently lagging orientation/area target. Update a pure
                // reaction face from the prescribed member and load vectors.
                V net=V::Zero();int reactionFace=-1;
                for(int f=0;f<int(c.faces.size());++f)for(auto source:c.faces[f].sources){
                    if(source.type==zForcePolyhedraFaceType::Reaction){if(c.faces[f].sources.size()==1)reactionFace=f;}
                    else net+=c.normals[f]*source.magnitude;
                }
                if(reactionFace>=0&&net.norm()>1e-12){auto& f=c.faces[reactionFace];c.normals[reactionFace]=-net.normalized();f.normal=publicVector(c.normals[reactionFace]);f.targetArea=net.norm()*s.settings.forceScale;f.sources[0].magnitude=net.norm();}
                c.report.converged=false;pursue(c,s.settings,s.settings.areaIterationsPerStep);
            }
        }
        s.moveForm(directions);++s.report.equilibriumIterations;s.summarize();
    }
    s.updateGraph();return s.report;
}
zForcePolyhedraReport zTsForcePolyhedra::solveEquilibrium(int maxIterations){
    require(maxIterations>0&&maxIterations<=100000,"Invalid equilibrium iteration count.");if(impl->cells.empty())createForcePolyhedra();return stepEquilibrium(maxIterations);
}
void zTsForcePolyhedra::reset(){impl->positions=impl->original;impl->precise.clear();for(auto p:impl->positions)impl->precise.push_back(vector(p));impl->invalidate();impl->updateGraph();}
const zObjectGraph& zTsForcePolyhedra::resultGraph()const{return impl->graph;}
const zIntArray& zTsForcePolyhedra::edgeConnects()const{return impl->edges;}
const zDoubleArray& zTsForcePolyhedra::axialForces()const{return impl->forces;}
const zVectorArray& zTsForcePolyhedra::residualForces()const{return impl->residuals;}
const zVectorArray& zTsForcePolyhedra::supportReactions()const{return impl->reactions;}
const zObjectMesh& zTsForcePolyhedra::forceCell(int id)const{return impl->cells[impl->checkedCell(id)].mesh;}
const std::vector<zForcePolyhedraFace>& zTsForcePolyhedra::cellFaces(int id)const{return impl->cells[impl->checkedCell(id)].faces;}
const zForcePolyhedraCellReport& zTsForcePolyhedra::cellReport(int id)const{return impl->cells[impl->checkedCell(id)].report;}
const zForcePolyhedraReport& zTsForcePolyhedra::report()const{return impl->report;}
}
