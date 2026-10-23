#include <zspace/zToolsets/fabrication/zTsParametricStairs.h>
#include <Eigen/Dense>
#include "zTsStairsLayout.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>
#include <limits>

namespace zSpace {
void buildStairsSurface(zStairsLayoutData&,const zPointArray&,const zIntArray&,const zIntArray&,const zParametricStairsSettings&,bool);
namespace {
using V = Eigen::Vector2d;
struct Arm { int neighbor, left=-1, right=-1; V direction; double angle; };
void require(bool value, const char* message) { if(!value) throw std::invalid_argument(message); }
}
void zTsParametricStairs::update(zObjectGraph& graph, const zParametricStairsSettings& s) {
    require(std::isfinite(s.width)&&s.width>0&&s.width<=1e6,"Stairs width must be positive and at most 1e6.");
    require(s.divisions==-1||(s.divisions>=2&&s.divisions<=256),"Legacy stairs divisions must be 2..256.");
    require(s.numDivs>=0&&s.numDivs<=5,"Stairs numDivs must be 0..5.");
    require(s.numSmoothDivs>=0&&s.numSmoothDivs<=5,"Stairs numSmoothDivs must be 0..5.");
    const bool legacy=s.divisions!=-1;
    const int longitudinal=legacy?s.divisions:1;
    require(std::isfinite(s.miterLimit)&&s.miterLimit>=1&&s.miterLimit<=100,"Stairs miter limit must be 1..100.");
    require(std::isfinite(s.vertexMass)&&s.vertexMass>0,"Stairs mass must be positive.");
    zPointArray input; zIntArray edges; zFnGraph fn(graph);fn.getVertexPositions(input);fn.getEdgeData(edges);
    require(input.size()>=2&&input.size()<=4096&&edges.size()>=2&&edges.size()<=8192,"Stairs graph size outside limits.");
    require(s.vertexWidths.empty()||s.vertexWidths.size()==input.size(),"Stairs widths must match graph vertex count.");
    auto width=[&](int v){double w=s.vertexWidths.empty()?s.width:s.vertexWidths[v];require(std::isfinite(w)&&w>0&&w<=1e6,"Invalid graph vertex width.");return w;};
    std::vector<std::vector<Arm>> arms(input.size());std::set<std::pair<int,int>> unique;
    for(auto p:input)require(std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z)&&std::abs(p.x)<=1e8&&std::abs(p.y)<=1e8&&std::abs(p.z)<=1e8,"Invalid stairs graph coordinate.");
    for(size_t e=0;e<edges.size();e+=2){int a=edges[e],b=edges[e+1];
        require(a>=0&&b>=0&&a<int(input.size())&&b<int(input.size())&&a!=b,"Invalid stairs graph edge.");
        require(unique.emplace(std::min(a,b),std::max(a,b)).second,"Duplicate stairs edge.");
        V d(input[b].x-input[a].x,input[b].y-input[a].y);require(d.norm()>1e-6,"Stairs edges need nonzero XY length.");d.normalize();
        arms[a].push_back({b,-1,-1,d,std::atan2(d.y(),d.x())});arms[b].push_back({a,-1,-1,-d,std::atan2(-d.y(),-d.x())});
    }
    // Each connected component needs end supports; isolated nodes are invalid.
    std::set<int> visited;
    for(int v=0;v<int(input.size());++v){require(!arms[v].empty(),"Stairs graph contains an isolated vertex.");if(visited.count(v))continue;
        std::vector<int> queue={v};visited.insert(v);int ends=0;
        while(!queue.empty()){int i=queue.back();queue.pop_back();if(arms[i].size()==1)++ends;for(auto a:arms[i])if(visited.insert(a.neighbor).second)queue.push_back(a.neighbor);}
        require(ends>=2,"Each stairs component needs at least two valence-one supports.");
    }
    std::shared_ptr<zStairsLayoutData> layout;
    if(s.generateSteps){
        layout=buildStairsLayout(input,edges,s);input=layout->graphPoints;
        double maxRise=0;for(double r:layout->actualRisers)maxRise=std::max(maxRise,r);
        for(size_t f=0;f<layout->flightPaths.size();++f){
            double h=std::abs(layout->flightFloors[f][1]-layout->flightFloors[f][0]);
            double len=layout->startReserves[f]+layout->usableRuns[f]+layout->endReserves[f];
            double run=std::max(0,layout->riserCounts[f]-1)*(s.flightTreadDepths.empty()?s.treadDepth:s.flightTreadDepths[f]);
            double start=layout->startReserves[f]+(layout->usableRuns[f]-run)*.5;
            if(len>1e-6)maxRise=std::max(maxRise,h*(start/len));
        }
        double offset=s.clearance+.5*s.structureThickness+maxRise;
        for(auto& p:input)p.z-=float(offset);
    }
    zPointArray positions;zIntArray sources,centers(input.size()),supports,counts,connects,faceGraphEdges;
    auto add=[&](zPoint p,int source){int id=int(positions.size());positions.push_back(p);sources.push_back(source);return id;};
    for(int v=0;v<int(input.size());++v){auto& list=arms[v];std::sort(list.begin(),list.end(),[](const Arm& a,const Arm& b){return a.angle<b.angle;});
        centers[v]=add(input[v],v);double h=width(v)*.5;
        auto corner=[&](V offset){auto p=input[v];p.x+=float(offset.x());p.y+=float(offset.y());return add(p,v);};
        if(list.size()==1){V n(-list[0].direction.y(),list[0].direction.x());list[0].left=corner(n*h);list[0].right=corner(-n*h);supports.insert(supports.end(),{list[0].left,centers[v],list[0].right});}
        else for(size_t i=0;i<list.size();++i){auto& a=list[i];auto& b=list[(i+1)%list.size()];
            double gap=b.angle-a.angle;if(gap<=0)gap+=2*3.141592653589793;
            require(gap>1e-5&&gap<2*3.141592653589793-1e-5,"Coincident stairs junction directions.");
            // Angular bisector of adjacent offset lines, including reflex wedges.
            double mid=a.angle+gap*.5;
            // Intersection: along the wedge bisector at half-width / sin(gap/2).
            V offset=V(std::cos(mid),std::sin(mid))*(h/std::sin(gap*.5));
            require(offset.norm()<=h*s.miterLimit,"Stairs junction too sharp; reduce width or edit the graph.");
            int id=corner(offset);a.left=id;b.right=id;
        }
    }
    auto arm=[&](int a,int b)->Arm&{for(auto& x:arms[a])if(x.neighbor==b)return x;throw std::logic_error("Missing stairs arm.");};
    auto buildRibbon=[&](int sections){
    require(sections<=10000,"Stair construction exceeds 10000 stations per graph edge.");
    for(size_t e=0;e<edges.size();e+=2){int a=edges[e],b=edges[e+1];auto& aa=arm(a,b);auto& bb=arm(b,a);
        int start[3]={aa.right,centers[a],aa.left},end[3]={bb.left,centers[b],bb.right};
        std::vector<std::array<int,3>> rows={{start[0],start[1],start[2]}};
        for(int i=1;i<sections;++i){double t=double(i)/sections;std::array<int,3> row;
            for(int k=0;k<3;++k){auto p=positions[start[k]],q=positions[end[k]];row[k]=add(zPoint(float(p.x*(1-t)+q.x*t),float(p.y*(1-t)+q.y*t),float(p.z*(1-t)+q.z*t)),-1);}rows.push_back(row);}
        rows.push_back({end[0],end[1],end[2]});
        for(int i=0;i<sections;++i)for(int k=0;k<2;++k){int ids[4]={rows[i][k],rows[i+1][k],rows[i+1][k+1],rows[i][k+1]};
            // Positive XY area at every corner rejects strips folded by overly wide junctions.
            for(int j=0;j<4;++j){auto p=positions[ids[j]],q=positions[ids[(j+1)%4]],r=positions[ids[(j+2)%4]];
                double cross=double(q.x-p.x)*(r.y-q.y)-double(q.y-p.y)*(r.x-q.x);
                require(cross>1e-10,"Stairs ribbon folds or collapses; reduce width or junction sharpness.");}
            counts.push_back(4);faceGraphEdges.push_back(int(e/2));connects.insert(connects.end(),ids,ids+4);
        }
    }
    };
    buildRibbon(longitudinal);
    zObjectMesh mesh;zFnMesh(mesh).create(positions,counts,connects);
    zObjectMesh lowPoly=mesh;
    if(!legacy&&s.numDivs>0){
        // Only insert transverse stations along each graph arm. The existing
        // right/centre/left rows keep the same two faces across its width.
        counts.clear();connects.clear();faceGraphEdges.clear();buildRibbon(1<<s.numDivs);
        zFnMesh(mesh).create(positions,counts,connects);
    }
    // End-cap lines stay fixed through refinement. Their new boundary vertices
    // inherit the same graph-node support rather than relying on vertex ordering.
    struct Cap { int graphVertex; zPoint a,b; };
    std::vector<Cap> caps;
    for(int v=0;v<int(input.size());++v)if(arms[v].size()==1)caps.push_back({v,positions[arms[v][0].left],positions[arms[v][0].right]});
    auto refine=[&](){
        zFnMesh fn(mesh);zPointArray before;fn.getVertexPositions(before);
        require(size_t(fn.numVertices())+fn.numEdges()+fn.numPolygons()<size_t(std::numeric_limits<int>::max()/8),"Stairs refinement exceeds representable solver indices.");
        if(layout){zIntArray c,f,child;fn.getPolygonData(f,c);for(size_t i=0;i<c.size();++i)for(int k=0;k<c[i];++k)child.push_back(faceGraphEdges[i]);faceGraphEdges=std::move(child);}
        fn.smoothMesh(1,false);
        zPointArray refined;fn.getVertexPositions(refined);
        require(refined.size()>=before.size(),"Stairs refinement lost source vertices.");
        // Supported original vertices are exact anchors, including cap midpoints.
        for(int id:supports)refined[id]=before[id];
        sources.resize(refined.size(),-1);std::set<int> fixed(supports.begin(),supports.end());
        for(int id=int(before.size());id<int(refined.size());++id)for(const auto& cap:caps){
            Eigen::Vector3d a(cap.a.x,cap.a.y,cap.a.z),b(cap.b.x,cap.b.y,cap.b.z),p(refined[id].x,refined[id].y,refined[id].z);
            auto d=b-a;double t=(p-a).dot(d)/d.squaredNorm();
            if(t>=-1e-6&&t<=1+1e-6&&(p-a-t*d).norm()<=1e-6*std::max(1.0,d.norm())){sources[id]=cap.graphVertex;fixed.insert(id);break;}
        }
        supports.assign(fixed.begin(),fixed.end());fn.setVertexPositions(refined);
    };
    zObjectMesh subdivided=mesh;
    if(!legacy)for(int i=0;i<s.numSmoothDivs;++i)refine();
    if(layout){for(int e:faceGraphEdges)layout->meshFaceFlights.push_back(layout->edgeFlights[e]);zFnMesh fn(mesh);zPointArray p;fn.getVertexPositions(p);zIntArray c,f;fn.getPolygonData(f,c);buildStairsSurface(*layout,p,c,f,s,true);lowerStairsBelowEnvelope(p,c,f,supports,*layout,.5*s.structureThickness,s.clearance);fn.setVertexPositions(p);layout->minimumClearance=stairsClearance(p,c,f,*layout,.5*s.structureThickness);}
    auto vault=std::make_unique<zTsVault>();vault->setInputMesh(mesh);vault->setSupports(supports);vault->setVertexMass(s.vertexMass);vault->setBestFitSettings(s.bestFit);vault->computeBestFitForceDensities();
    m_mesh=std::move(mesh);m_lowPoly=std::move(lowPoly);m_subdivided=std::move(subdivided);m_sources=std::move(sources);m_centers=std::move(centers);m_supports=std::move(supports);m_vault=std::move(vault);m_settings=s;m_layout=std::move(layout);
    if(m_layout){m_designPoints=m_layout->graphPoints;m_designEdges=edges;m_flightCounts=m_layout->riserCounts;zFnMesh(m_mesh).getVertexPositions(m_centerReferences);attachClearanceValidator();}
}
zVaultSolveReport zTsParametricStairs::stepBestFit(int iterations){if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");auto r=m_vault->stepBestFit(iterations);refreshLayout();return r;}
zVaultSolveReport zTsParametricStairs::snapFixedDensities(double geometryWeight){if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");auto r=m_vault->snapFixedDensities(geometryWeight);refreshLayout();return r;}
void zTsParametricStairs::smoothResult(int divisions) {
    require(bool(m_vault),"Stairs graph is not initialized.");
    require(divisions>=1&&divisions<=5,"Post-FDM smooth divisions must be 1..5.");
    zObjectMesh mesh=resultMesh();zIntArray supports=m_supports,faceFlights=m_layout?m_layout->meshFaceFlights:zIntArray{};
    for(int pass=0;pass<divisions;++pass){
        zFnMesh fn(mesh);zPointArray before;fn.getVertexPositions(before);
        std::set<int> fixed(supports.begin(),supports.end());
        std::vector<std::pair<zPoint,zPoint>> caps;
        for(zItMeshEdge e(mesh);!e.end();e++) {zIntArray ids;e.getVertices(ids);if(fixed.count(ids[0])&&fixed.count(ids[1]))caps.emplace_back(before[ids[0]],before[ids[1]]);}
        require(size_t(fn.numVertices())+fn.numEdges()+fn.numPolygons()<size_t(std::numeric_limits<int>::max()/8),"Smoothing exceeds representable indices.");
        if(m_layout){zIntArray c,f,child;fn.getPolygonData(f,c);for(size_t i=0;i<c.size();++i)for(int k=0;k<c[i];++k)child.push_back(faceFlights[i]);faceFlights=std::move(child);}
        fn.smoothMesh(1,false);zPointArray after;fn.getVertexPositions(after);
        for(int v:supports)after[v]=before[v];
        // Refinement must not shorten flights whose exact treads already fill the run.
        // Keep editing-node plan anchors; free heights still receive smoothing.
        if(m_layout)for(int v:m_centers){after[v].x=before[v].x;after[v].y=before[v].y;}
        for(int v=int(before.size());v<int(after.size());++v)for(const auto& cap:caps){
            Eigen::Vector3d a(cap.first.x,cap.first.y,cap.first.z),b(cap.second.x,cap.second.y,cap.second.z),p(after[v].x,after[v].y,after[v].z);
            auto d=b-a;double length=d.squaredNorm();if(length<=1e-20)continue;
            double t=(p-a).dot(d)/length;
            if(t>=-1e-6&&t<=1+1e-6&&(p-a-t*d).norm()<=1e-6*std::max(1.0,d.norm())){fixed.insert(v);break;}
        }
        supports.assign(fixed.begin(),fixed.end());fn.setVertexPositions(after);
    }
    std::shared_ptr<zStairsLayoutData> layout;
    if(m_layout){zFnMesh fn(mesh);zPointArray p;fn.getVertexPositions(p);layout=layoutFor(p,true);layout->meshFaceFlights=faceFlights;zIntArray c,f;fn.getPolygonData(f,c);buildStairsSurface(*layout,p,c,f,m_settings,true);lowerStairsBelowEnvelope(p,c,f,supports,*layout,.5*m_settings.structureThickness,m_settings.clearance);fn.setVertexPositions(p);layout->minimumClearance=stairsClearance(p,c,f,*layout,.5*m_settings.structureThickness);}
    auto vault=std::make_unique<zTsVault>();vault->setInputMesh(mesh);vault->setSupports(supports);vault->setVertexMass(m_settings.vertexMass);vault->setBestFitSettings(m_settings.bestFit);vault->computeBestFitForceDensities();
    m_mesh=std::move(mesh);m_supports=std::move(supports);m_sources.resize(zFnMesh(m_mesh).numVertices(),-1);m_vault=std::move(vault);if(layout){m_layout=std::move(layout);attachClearanceValidator();}
}
void zTsParametricStairs::resetBestFit(){if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");m_vault->resetBestFit();refreshLayout();}
bool zTsParametricStairs::previousBestFit(){if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");bool r=m_vault->previousBestFit();refreshLayout();return r;}
zObjectMesh& zTsParametricStairs::resultMesh(){if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->resultMesh();}
std::shared_ptr<zStairsLayoutData> zTsParametricStairs::layoutFor(const zPointArray& p,bool geometry) const {
    zPointArray graph=m_designPoints;
    for(size_t v=0;v<graph.size();++v){int i=m_centers[v];graph[v].x+=p[i].x-m_centerReferences[i].x;graph[v].y+=p[i].y-m_centerReferences[i].y;}
    auto contains=[](const zPoint& point,const zPointArray& poly){bool in=false;Eigen::Vector2d p(point.x,point.y);for(size_t i=0,j=poly.size()-1;i<poly.size();j=i++){auto a=poly[j],b=poly[i];Eigen::Vector2d pa(a.x,a.y),pb(b.x,a.y);pb.y()=b.y;auto d=pb-pa;double l2=d.squaredNorm();if(l2>1e-20&&(p-pa-std::clamp((p-pa).dot(d)/l2,0.0,1.0)*d).squaredNorm()<=1e-6)return true;if((a.y>point.y)!=(b.y>point.y)&&point.x<(b.x-a.x)*(point.y-a.y)/(b.y-a.y)+a.x)in=!in;}return in;};
    for(const auto& poly:m_settings.landingPolygons)for(size_t v=0;v<graph.size();++v)if(std::abs(m_designPoints[v].z-poly.front().z)<1e-3&&contains(m_designPoints[v],poly))require(contains(graph[v],poly),"A solved graph node left its fixed Rhino landing footprint.");
    auto layout=buildStairsLayout(graph,m_designEdges,m_settings,m_flightCounts,geometry);layout->meshFaceFlights=m_layout->meshFaceFlights;return layout;
}
void zTsParametricStairs::attachClearanceValidator(){
    zIntArray c,f;zFnMesh(resultMesh()).getPolygonData(f,c);
    m_vault->setGeometryValidator([this,c,f](const zPointArray& p){try{auto l=layoutFor(p);buildStairsSurface(*l,p,c,f,m_settings,false);return stairsClearance(p,c,f,*l,.5*m_settings.structureThickness)>=m_settings.clearance-1e-6;}catch(const std::exception&){return false;}});
}
void zTsParametricStairs::refreshLayout(){
    if(!m_layout)return;zFnMesh fn(resultMesh());zPointArray p;fn.getVertexPositions(p);auto layout=layoutFor(p,true);zIntArray c,f;fn.getPolygonData(f,c);buildStairsSurface(*layout,p,c,f,m_settings,true);layout->minimumClearance=stairsClearance(p,c,f,*layout,.5*m_settings.structureThickness);m_layout=std::move(layout);
}
zObjectMesh& zTsParametricStairs::stairMesh(){if(!m_layout)throw std::logic_error("Generate stair layout first.");return m_layout->mesh;}
zObjectGraph& zTsParametricStairs::walkingGraph(){if(!m_layout)throw std::logic_error("Generate stair layout first.");return m_layout->graph;}
const zIntArray& zTsParametricStairs::stairFaceTypes()const{if(!m_layout)throw std::logic_error("Generate stair layout first.");return m_layout->faceTypes;}
const zIntArray& zTsParametricStairs::flightRiserCounts()const{if(!m_layout)throw std::logic_error("Generate stair layout first.");return m_layout->riserCounts;}
const zDoubleArray& zTsParametricStairs::flightActualRisers()const{if(!m_layout)throw std::logic_error("Generate stair layout first.");return m_layout->actualRisers;}
double zTsParametricStairs::minimumClearance()const{return m_layout?m_layout->minimumClearance:0;}
const zDoubleArray& zTsParametricStairs::forceDensities() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->forceDensities();}
const zIntArray& zTsParametricStairs::edgeConnects() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->edgeConnects();}
const zIntArray& zTsParametricStairs::activeEdgeIds() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->activeEdgeIds();}
const zDoubleArray& zTsParametricStairs::axialForces() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->axialForces();}
const zVectorArray& zTsParametricStairs::residualForces() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->residualForces();}
const zVectorArray& zTsParametricStairs::supportReactions() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->supportReactions();}
const zVaultSolveReport& zTsParametricStairs::report() const {if(!m_vault)throw std::logic_error("Stairs graph is not initialized.");return m_vault->report();}
}
