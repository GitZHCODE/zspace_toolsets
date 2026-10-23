#include "zTsStairsLayout.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <stdexcept>
#include <limits>

namespace zSpace {
namespace {
using V=Eigen::Vector2d;
V xy(const zPoint& p){return {p.x,p.y};}
double cross(V a,V b){return a.x()*b.y()-a.y()*b.x();}
void require(bool b,const char* m){if(!b)throw std::invalid_argument(m);}
bool inside(V p,const zPointArray& poly,double tol=1e-7){
    bool in=false;for(size_t i=0,j=poly.size()-1;i<poly.size();j=i++){
        V a=xy(poly[j]),b=xy(poly[i]),d=b-a;
        if(d.squaredNorm()>1e-20&&(p-a-std::clamp((p-a).dot(d)/d.squaredNorm(),0.0,1.0)*d).squaredNorm()<=tol*tol)return true;
        if((a.y()>p.y())!=(b.y()>p.y())&&p.x()<(b.x()-a.x())*(p.y()-a.y())/(b.y()-a.y())+a.x())in=!in;
    }return in;
}
zPointArray hull(zPointArray p){
    std::sort(p.begin(),p.end(),[](auto a,auto b){return a.x<b.x||(a.x==b.x&&a.y<b.y);});
    zPointArray out;for(auto a:p){while(out.size()>1&&cross(xy(out.back())-xy(out[out.size()-2]),xy(a)-xy(out.back()))<=1e-9)out.pop_back();out.push_back(a);}
    size_t n=out.size();for(int i=int(p.size())-2;i>=0;--i){while(out.size()>n&&cross(xy(out.back())-xy(out[out.size()-2]),xy(p[i])-xy(out.back()))<=1e-9)out.pop_back();out.push_back(p[i]);}
    if(out.size()>1)out.pop_back();return out;
}
std::vector<V> overlap(const std::array<zPoint,3>& tri,const zPointArray& polygon){
    std::vector<V> out;const V a=xy(tri[0]),b=xy(tri[1]),c=xy(tri[2]);const double area=cross(b-a,c-a);
    auto inTri=[&](V p){double u=cross(b-p,c-p)/area,v=cross(c-p,a-p)/area,w=1-u-v;return u>=-1e-7&&v>=-1e-7&&w>=-1e-7;};
    for(auto p:tri)if(inside(xy(p),polygon))out.push_back(xy(p));
    for(auto p:polygon)if(inTri(xy(p)))out.push_back(xy(p));
    for(int i=0;i<3;++i)for(size_t j=0;j<polygon.size();++j){V p=xy(tri[i]),d=xy(tri[(i+1)%3])-p,q=xy(polygon[j]),e=xy(polygon[(j+1)%polygon.size()])-q;double denom=cross(d,e);if(std::abs(denom)<1e-12)continue;
        double t=cross(q-p,e)/denom,u=cross(q-p,d)/denom;if(t>=-1e-7&&t<=1+1e-7&&u>=-1e-7&&u<=1+1e-7)out.push_back(p+std::clamp(t,0.0,1.0)*d);
    }return out;
}
struct Bounds {double x0=1e100,y0=1e100,x1=-1e100,y1=-1e100;};
Bounds bounds(const zPointArray& p){Bounds b;for(auto v:p){b.x0=std::min(b.x0,double(v.x));b.x1=std::max(b.x1,double(v.x));b.y0=std::min(b.y0,double(v.y));b.y1=std::max(b.y1,double(v.y));}return b;}
bool disjoint(Bounds a,Bounds b){return a.x1<b.x0-1e-7||b.x1<a.x0-1e-7||a.y1<b.y0-1e-7||b.y1<a.y0-1e-7;}
template<class F> void intersections(const zPointArray& points,const zIntArray& counts,const zIntArray& connects,const zStairsLayoutData& layout,F visit){
    std::vector<Bounds> boxes;for(auto& f:layout.footprints)boxes.push_back(bounds(f.polygon));
    // Index projected envelopes so refinement does not require all-pairs tests.
    Bounds all;for(auto b:boxes){all.x0=std::min(all.x0,b.x0);all.y0=std::min(all.y0,b.y0);all.x1=std::max(all.x1,b.x1);all.y1=std::max(all.y1,b.y1);}
    double cell=std::max(1e-5,std::max(all.x1-all.x0,all.y1-all.y0)/std::max(1.0,std::sqrt(double(boxes.size()))));
    auto range=[&](Bounds b){return std::array<int,4>{int(std::floor((b.x0-all.x0)/cell)),int(std::floor((b.y0-all.y0)/cell)),int(std::floor((b.x1-all.x0)/cell)),int(std::floor((b.y1-all.y0)/cell))};};
    std::map<std::pair<int,int>,std::vector<size_t>> grid;std::vector<size_t> large;
    for(size_t f=0;f<boxes.size();++f){auto r=range(boxes[f]);if(double(r[2]-r[0]+1)*(r[3]-r[1]+1)>4096){large.push_back(f);continue;}for(int x=r[0];x<=r[2];++x)for(int y=r[1];y<=r[3];++y)grid[{x,y}].push_back(f);}
    std::vector<size_t> seen(boxes.size(),0);size_t query=0;
    size_t start=0,face=0;for(int count:counts){int owner=layout.meshFaceFlights.empty()?-1:layout.meshFaceFlights.at(face);for(int i=1;i+1<count;++i){std::array<int,3> ids={connects[start],connects[start+i],connects[start+i+1]};std::array<zPoint,3> tri={points[ids[0]],points[ids[1]],points[ids[2]]};
        double area=cross(xy(tri[1])-xy(tri[0]),xy(tri[2])-xy(tri[0]));require(area>1e-10,"Stairs structural mesh has a folded or collapsed plan triangle.");Bounds box=bounds({tri[0],tri[1],tri[2]});
        auto r=range(box);std::vector<size_t> candidates=large;++query;
        if(double(r[2]-r[0]+1)*(r[3]-r[1]+1)>4096){candidates.clear();for(size_t f=0;f<boxes.size();++f)candidates.push_back(f);}else for(int x=r[0];x<=r[2];++x)for(int y=r[1];y<=r[3];++y){auto it=grid.find({x,y});if(it!=grid.end())candidates.insert(candidates.end(),it->second.begin(),it->second.end());}
        for(size_t f:candidates){if(seen[f]==query)continue;seen[f]=query;if(disjoint(box,boxes[f]))continue;std::vector<std::array<double,3>> weights;for(V p:overlap(tri,layout.footprints[f].polygon)){
            double u=cross(xy(tri[1])-p,xy(tri[2])-p)/area,v=cross(xy(tri[2])-p,xy(tri[0])-p)/area;weights.push_back({u,v,1-u-v});
        }if(!weights.empty())visit(ids,weights,layout.footprints[f],owner);}
    }start+=count;++face;}
}
bool ownedFootprint(const zStairsLayoutData& layout,const zStairsFootprint& footprint,int owner){
    if(owner<0)return true;if(footprint.flight==owner)return true;
    if(footprint.flight<0){auto floors=layout.flightFloors.at(owner);return std::abs(footprint.top-floors[0])<1e-4||std::abs(footprint.top-floors[1])<1e-4;}
    return false;
}
}
std::shared_ptr<zStairsLayoutData> buildStairsLayout(const zPointArray& points,const zIntArray& edges,const zParametricStairsSettings& s,const zIntArray& lockedCounts,bool geometry,bool capacityOnly){
    for(double v:{s.treadDepth,s.riserHeight,s.landingDepth,s.structureThickness})require(std::isfinite(v)&&v>0,"Stair tread, riser, landing and thickness inputs must be positive.");
    require(std::isfinite(s.clearance)&&s.clearance>=0,"Stair clearance must be nonnegative.");
    auto out=std::make_shared<zStairsLayoutData>();out->graphPoints=points;
    std::vector<std::vector<int>> neighbors(points.size());
    for(size_t i=0;i<edges.size();i+=2){int a=edges[i],b=edges[i+1];require(a>=0&&b>=0&&a<int(points.size())&&b<int(points.size())&&a!=b,"Invalid stairs graph edge.");neighbors[a].push_back(b);neighbors[b].push_back(a);}
    std::vector<zPointArray> landings(points.size());out->edgeFlights.assign(edges.size()/2,-1);
    auto width=[&](int v){return s.vertexWidths.empty()?s.width:s.vertexWidths.at(v);};
    for(auto poly:s.landingPolygons){if(poly.size()>2&&(poly.front()-poly.back()).length()<1e-6)poly.pop_back();require(poly.size()>=3,"Landing outline needs at least three points.");double z=poly.front().z;for(auto p:poly)require(std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z)&&std::abs(p.z-z)<1e-4,"Landing outlines must be horizontal.");
        for(int v=0;v<int(points.size());++v)if(inside(xy(points[v]),poly,1e-3)){if(std::abs(points[v].z-z)>=1e-3)continue;require(landings[v].empty(),"Multiple landing outlines contain the same graph node.");landings[v]=poly;}
    }
    for(int v:s.landingVertices)require(v>=0&&v<int(points.size()),"Invalid selected landing vertex.");
    auto node=[&](int v){return neighbors[v].size()!=2||(s.landingVertices.empty()?!landings[v].empty():std::find(s.landingVertices.begin(),s.landingVertices.end(),v)!=s.landingVertices.end());};
    for(int v=0;v<int(points.size());++v)if(node(v)&&landings[v].empty()){
        zPointArray corners;for(int n:neighbors[v]){V d=xy(points[n])-xy(points[v]);require(d.norm()>1e-6,"Collapsed stair flight plan.");d.normalize();V normal(-d.y(),d.x());
            for(double along:{-.5*s.landingDepth,.5*s.landingDepth})for(double side:{-.5*width(v),.5*width(v)}){V p=xy(points[v])+along*d+side*normal;corners.emplace_back(float(p.x()),float(p.y()),points[v].z);}}
        landings[v]=hull(corners);
    }
    int currentFlight=-1;
    zPointArray vertices;zIntArray counts,connects;
    auto face=[&](zPointArray p,int type){if(!geometry)return;counts.push_back(int(p.size()));for(auto v:p){connects.push_back(int(vertices.size()));vertices.push_back(v);}out->faceTypes.push_back(type);};
    auto slab=[&](const zPointArray& poly,double top,int type){zPointArray upper=poly;for(auto& p:upper)p.z=float(top);
        out->footprints.push_back({upper,top,top,currentFlight});face(upper,type);
    };
    // An imported polygon may cover several nodes; publish each footprint once.
    std::set<std::vector<float>> published;
    for(auto poly:landings)if(!poly.empty()){std::vector<float> key;for(auto p:poly)key.insert(key.end(),{p.x,p.y,p.z});if(published.insert(key).second)slab(poly,poly.front().z,2);}
    for(auto poly:s.landingPolygons){if(poly.size()>2&&(poly.front()-poly.back()).length()<1e-6)poly.pop_back();std::vector<float> key;for(auto p:poly)key.insert(key.end(),{p.x,p.y,p.z});if(published.insert(key).second)slab(poly,poly.front().z,2);}
    std::set<std::pair<int,int>> used;size_t flight=0;
    for(int v=0;v<int(points.size());++v)if(node(v))for(int next:neighbors[v]){
        if(used.count(std::minmax(v,next)))continue;std::vector<int> path={v};int previous=v,current=next;
        while(true){require(used.emplace(std::min(previous,current),std::max(previous,current)).second,"Stair flight loop is unsupported.");path.push_back(current);if(node(current))break;int following=neighbors[current][0]==previous?neighbors[current][1]:neighbors[current][0];previous=current;current=following;require(path.size()<=points.size()+1,"Stair flight loop is unsupported.");}
        require(path.front()!=path.back(),"A stair flight cannot return to the same landing.");
        if(points[path.front()].z>points[path.back()].z)std::reverse(path.begin(),path.end());
        std::vector<double> station={0};for(size_t i=1;i<path.size();++i)station.push_back(station.back()+(xy(points[path[i]])-xy(points[path[i-1]])).norm());double length=station.back();
        auto reserve=[&](int a,int b){V origin=xy(points[a]),d=xy(points[b])-origin;d.normalize();double exit=std::numeric_limits<double>::infinity();auto& poly=landings[a];
            for(size_t i=0;i<poly.size();++i){V p=xy(poly[i]),edge=xy(poly[(i+1)%poly.size()])-p;double denom=cross(d,edge);if(std::abs(denom)<1e-10)continue;double t=cross(p-origin,edge)/denom,u=cross(p-origin,d)/denom;if(t>=-1e-3&&u>=-1e-6&&u<=1+1e-6)exit=std::min(exit,std::max(0.0,t));}
            require(std::isfinite(exit),"Graph node lies outside its landing outline.");return exit;
        };
        double outlineA=reserve(path.front(),path[1]),outlineB=reserve(path.back(),path[path.size()-2]);
        double a=std::max(.5*s.landingDepth,outlineA),b=std::max(.5*s.landingDepth,outlineB),height=double(points[path.back()].z)-points[path.front()].z;
        double tread=s.flightTreadDepths.empty()?s.treadDepth:s.flightTreadDepths.at(flight);
        require(std::isfinite(tread)&&tread>0,"Invalid selected flight tread depth.");
        double desired=height<=1e-7?0:std::ceil(height/s.riserHeight-1e-9);
        double capacity=std::floor(std::max(0.0,length-a-b)/tread+1e-7)+1;
        require((capacityOnly?capacity:desired)<=10000,"Stair flight exceeds 10000 risers.");
        int n=capacityOnly?int(capacity):!lockedCounts.empty()?lockedCounts.at(flight):!s.flightRiserCounts.empty()?s.flightRiserCounts.at(flight):int(desired);
        require(n>=0&&n<=10000,"Stair flight exceeds 10000 risers.");double rise=n?height/n:0,run=std::max(0,n-1)*tread;
        if(length-a-b<run-1e-5&&!capacityOnly)throw std::invalid_argument("Stair flight "+std::to_string(flight+1)+" is too short: "+std::to_string(length-a-b)+" usable run, "+std::to_string(run)+" required for "+std::to_string(n)+" risers. Edit the graph, tread depth or landing space.");
        out->flightPaths.push_back(zIntArray(path.begin(),path.end()));out->usableRuns.push_back(length-a-b);out->startReserves.push_back(a);out->endReserves.push_back(b);
        if(capacityOnly){out->riserCounts.push_back(n);out->actualRisers.push_back(rise);++flight;continue;}
        require((n==0)==(height<=1e-7),"Selected riser count must agree with whether the flight has a rise.");
        // Preserve tread depth: unused run is shared by flat approaches to the landings.
        double start=a+(length-a-b-run)*.5,end=start+run;currentFlight=int(flight);
        for(size_t i=1;i<path.size();++i)for(size_t e=0;e<edges.size();e+=2)if(std::minmax(path[i-1],path[i])==std::minmax(edges[e],edges[e+1]))out->edgeFlights[e/2]=currentFlight;
        out->flightFloors.push_back({points[path.front()].z,points[path.back()].z});
        out->riserCounts.push_back(n);out->actualRisers.push_back(rise);++flight;
        auto sample=[&](double t){size_t i=1;while(i+1<station.size()&&station[i]<t)i++;double f=(t-station[i-1])/(station[i]-station[i-1]);auto p=points[path[i-1]],q=points[path[i]];return zPoint(float(p.x+(q.x-p.x)*f),float(p.y+(q.y-p.y)*f),0);};
        auto sectionOffset=[&](double t){
            size_t i=1;while(i+1<station.size()&&station[i]<t-1e-7)i++;
            double f=(t-station[i-1])/(station[i]-station[i-1]);double h=.5*(width(path[i-1])*(1-f)+width(path[i])*f);
            V d=xy(points[path[i]])-xy(points[path[i-1]]);d.normalize();V normal(-d.y(),d.x());
            if(i+1<station.size()&&std::abs(t-station[i])<1e-6){V next=xy(points[path[i+1]])-xy(points[path[i]]);next.normalize();V other(-next.y(),next.x());double denominator=1+normal.dot(other);require(denominator>1e-8,"Stair flight turn is too sharp.");V miter=(normal+other)*h/denominator;require(miter.norm()<=h*s.miterLimit,"Stair tread miter exceeds the junction limit.");return V(miter);}
            return V(normal*h);
        };
        auto strip=[&](double from,double to,double elevation,int type){if(to-from<1e-7)return;
            std::vector<double> cuts={from};for(double t:station)if(t>from+1e-7&&t<to-1e-7)cuts.push_back(t);cuts.push_back(to);
            for(size_t j=1;j<cuts.size();++j){auto p=sample(cuts[j-1]),q=sample(cuts[j]);V offset0=sectionOffset(cuts[j-1]),offset1=sectionOffset(cuts[j]);
                zPointArray poly;for(auto pair:{std::pair<zPoint,V>{p,-offset0},{q,-offset1},{q,offset1},{p,offset0}}){V v=xy(pair.first)+pair.second;poly.emplace_back(float(v.x()),float(v.y()),float(elevation));}
                for(size_t k=0;k<poly.size();++k)require(cross(xy(poly[(k+1)%4])-xy(poly[k]),xy(poly[(k+2)%4])-xy(poly[(k+1)%4]))>1e-10,"Stair tread folds around a flight turn; reduce width or edit the graph.");
                slab(poly,elevation,type);
            }
        };
        strip(outlineA,start,points[path.front()].z,2);strip(end,length-outlineB,points[path.back()].z,2);
        if(!n)strip(start,end,points[path.front()].z,2);
        for(int k=0;k<n-1;++k)strip(start+k*tread,start+(k+1)*tread,points[path.front()].z+(k+1)*rise,0);
        for(size_t i=1;i+1<path.size();++i){double t=station[i];double f=run>1e-9?std::clamp((t-start)/run,0.0,1.0):t>=start?1:0;out->graphPoints[path[i]].z=float(points[path.front()].z+height*f);}
    }
    require(!out->riserCounts.empty(),"Stair layout needs flights between landings.");
    require(s.flightTreadDepths.empty()||flight==s.flightTreadDepths.size(),"Selected flight tread count does not match the graph.");
    require(s.flightRiserCounts.empty()||flight==s.flightRiserCounts.size(),"Selected flight riser count does not match the graph.");
    if(!lockedCounts.empty())require(flight==lockedCounts.size(),"Stair flight topology changed during form finding.");
    if(geometry){zFnMesh(out->mesh).create(vertices,counts,connects);auto gp=out->graphPoints;auto ge=edges;zFnGraph(out->graph).create(gp,ge);}
    return out;
}
zStairsGraphOptimisation zTsParametricStairs::optimiseGraph(zObjectGraph& graph,
    const zParametricStairsSettings& input,double treadMin,double treadMax,double riserMin,double riserMax){
    require(std::isfinite(treadMin)&&std::isfinite(treadMax)&&treadMin>0&&treadMax>=treadMin,"Require 0 < tread minimum <= tread maximum.");
    require(std::isfinite(riserMin)&&std::isfinite(riserMax)&&riserMin>0&&riserMax>=riserMin,"Require 0 < riser minimum <= riser maximum.");
    zPointArray original;zIntArray edges;zFnGraph fn(graph);fn.getVertexPositions(original);fn.getEdgeData(edges);
    require(original.size()>=2&&original.size()<=4096&&edges.size()>=2&&edges.size()<=8192,"Stairs graph size outside limits.");
    require(input.vertexWidths.empty()||input.vertexWidths.size()==original.size(),"Stairs widths must match graph vertex count.");
    std::vector<zIntArray> adjacent(original.size());std::set<std::pair<int,int>> unique;
    for(size_t v=0;v<original.size();++v){auto p=original[v];require(std::isfinite(p.x)&&std::isfinite(p.y)&&std::isfinite(p.z),"Invalid graph coordinate.");double w=input.vertexWidths.empty()?input.width:input.vertexWidths[v];require(std::isfinite(w)&&w>0,"Invalid graph width.");}
    for(size_t e=0;e<edges.size();e+=2){int a=edges[e],b=edges[e+1];require(a>=0&&b>=0&&a<int(original.size())&&b<int(original.size())&&a!=b,"Invalid stairs graph edge.");require(unique.emplace(std::min(a,b),std::max(a,b)).second,"Duplicate stairs graph edge.");adjacent[a].push_back(b);adjacent[b].push_back(a);}
    zStairsGraphOptimisation result;auto settings=input;settings.flightRiserCounts.clear();settings.flightTreadDepths.clear();settings.treadDepth=treadMin;
    std::set<int> visited;
    for(int v=0;v<int(original.size());++v){require(!adjacent[v].empty(),"Stairs graph contains an isolated vertex.");if(visited.count(v))continue;zIntArray queue={v};visited.insert(v);int bottom=v;
        for(size_t i=0;i<queue.size();++i){int u=queue[i];if(original[u].z<original[bottom].z-1e-7||(std::abs(original[u].z-original[bottom].z)<=1e-7&&adjacent[bottom].size()==2&&adjacent[u].size()!=2))bottom=u;for(int n:adjacent[u])if(visited.insert(n).second)queue.push_back(n);}
        result.bottomNodes.push_back(bottom);
        bool landing=adjacent[bottom].size()!=2;for(const auto& poly:settings.landingPolygons)if(!poly.empty()&&std::abs(poly.front().z-original[bottom].z)<1e-3&&inside(xy(original[bottom]),poly,1e-3))landing=true;
        if(!landing){zPointArray corners;double width=settings.vertexWidths.empty()?settings.width:settings.vertexWidths[bottom];for(int n:adjacent[bottom]){V d=xy(original[n])-xy(original[bottom]);require(d.norm()>1e-6,"Collapsed graph plan edge.");d.normalize();V normal(-d.y(),d.x());for(double along:{-.5*settings.landingDepth,.5*settings.landingDepth})for(double side:{-.5*width,.5*width}){V p=xy(original[bottom])+along*d+side*normal;corners.emplace_back(float(p.x()),float(p.y()),original[bottom].z);}}settings.landingPolygons.push_back(hull(corners));}
    }
    auto capacity=buildStairsLayout(original,edges,settings,{},false,true);
    size_t count=capacity->flightPaths.size();std::vector<zIntArray> flightsAt(original.size());
    for(size_t f=0;f<count;++f){const auto& path=capacity->flightPaths[f];require(capacity->usableRuns[f]>=-1e-5,"Minimum landing depth leaves no room for a stair flight; reduce landing depth or edit the graph.");flightsAt[path.front()].push_back(int(f));flightsAt[path.back()].push_back(int(f));}
    for(int v=0;v<int(original.size());++v)if(!flightsAt[v].empty())result.landingVertices.push_back(v);
    result.riserCounts.assign(count,-1);result.treadDepths.resize(count);result.riserHeights.resize(count);
    auto points=original;std::vector<bool> assigned(points.size(),false);
    struct Choice {int n=-1;double tread=0,riser=0,level=0,distance=1e100,score=1e100;};
    const double targetTread=std::clamp(input.treadDepth,treadMin,treadMax),targetRiser=std::clamp(input.riserHeight,riserMin,riserMax);
    auto choose=[&](size_t f,double parent,double desired,double bottom,bool modify){Choice best;double usable=std::max(0.0,capacity->usableRuns[f]);int maxN=capacity->riserCounts[f];
        for(int n=1;n<=maxN;++n){double maxT=n==1?treadMax:std::min(treadMax,usable/(n-1));if(maxT<treadMin-1e-7)continue;double tread=std::clamp(targetTread,treadMin,std::max(treadMin,maxT));
            for(int sign:{1,-1}){double maxR=sign<0?std::min(riserMax,(parent-bottom)/n):riserMax;if(maxR<riserMin-1e-7)continue;double required=sign*(desired-parent)/n;
                if(!modify&&(required<riserMin-1e-7||required>maxR+1e-7))continue;
                double rise=modify?std::clamp(required,riserMin,std::max(riserMin,maxR)):required;double level=parent+sign*n*rise;
                double distance=modify?std::abs(level-desired):0;double score=std::pow((rise-targetRiser)/std::max(1e-6,riserMax-riserMin),2)+std::pow((tread-targetTread)/std::max(1e-6,treadMax-treadMin),2);
                if(distance<best.distance-1e-7||(std::abs(distance-best.distance)<1e-7&&(score<best.score-1e-7||(std::abs(score-best.score)<1e-7&&n>best.n))))best={n,tread,rise,level,distance,score};
            }
        }return best;
    };
    for(int bottomNode:result.bottomNodes){double bottom=original[bottomNode].z;zIntArray queue={bottomNode};assigned[bottomNode]=true;
        for(size_t i=0;i<queue.size();++i){int parent=queue[i];for(int f:flightsAt[parent]){if(result.riserCounts[f]>=0)continue;const auto& path=capacity->flightPaths[f];int child=path.front()==parent?path.back():path.front();bool flat=std::abs(original[path.front()].z-original[path.back()].z)<=1e-7;Choice choice;
                if(flat){choice={0,targetTread,0,points[parent].z,0,0};if(assigned[child])require(std::abs(points[child].z-points[parent].z)<=1e-5,"A graph loop requires conflicting landing levels; adjust the graph or dimension ranges.");}
                else {double desired=assigned[child]?points[child].z:original[child].z;choice=choose(f,points[parent].z,desired,bottom,false);if(choice.n<0&&!assigned[child])choice=choose(f,points[parent].z,desired,bottom,true);require(choice.n>=0,"A graph loop cannot fit the selected tread/riser ranges at consistent landing levels.");}
                result.riserCounts[f]=choice.n;result.treadDepths[f]=choice.tread;result.riserHeights[f]=choice.riser;
                if(!assigned[child]){points[child].z=float(choice.level);assigned[child]=true;queue.push_back(child);}
        }}
    }
    // Move imported landing outlines only with their associated nodes, including stacked levels.
    result.landingPolygons=settings.landingPolygons;result.landingOffsets.assign(settings.landingPolygons.size(),0);
    for(size_t i=0;i<settings.landingPolygons.size();++i){const auto& poly=settings.landingPolygons[i];bool found=false;double delta=0;
        for(size_t v=0;v<original.size();++v)if(!poly.empty()&&std::abs(original[v].z-poly.front().z)<1e-3&&inside(xy(original[v]),poly,1e-3)){double d=points[v].z-original[v].z;if(found)require(std::abs(d-delta)<1e-4,"One Rhino landing outline spans nodes needing different levels.");found=true;delta=d;}
        result.landingOffsets[i]=delta;for(auto& p:result.landingPolygons[i])p.z+=float(delta);
    }
    for(size_t f=0;f<count;++f){const auto& path=capacity->flightPaths[f];double height=double(points[path.back()].z)-points[path.front()].z,usable=capacity->usableRuns[f],run=std::max(0,result.riserCounts[f]-1)*result.treadDepths[f],start=capacity->startReserves[f]+(usable-run)*.5;double station=0;
        for(size_t i=1;i+1<path.size();++i){station+=(xy(points[path[i]])-xy(points[path[i-1]])).norm();double t=run>1e-9?std::clamp((station-start)/run,0.0,1.0):station>=start?1:0;points[path[i]].z=float(points[path.front()].z+height*t);}
        result.riserHeights[f]=result.riserCounts[f]?std::abs(height)/result.riserCounts[f]:0;
    }
    for(size_t v=0;v<points.size();++v)if(!flightsAt[v].empty()&&std::abs(points[v].z-original[v].z)>1e-5)++result.adjustedLandings;
    settings.landingPolygons=result.landingPolygons;settings.landingVertices=result.landingVertices;settings.flightRiserCounts=result.riserCounts;settings.flightTreadDepths=result.treadDepths;settings.treadDepth=input.treadDepth;
    buildStairsLayout(points,edges,settings,{},false); // Validate the whole selected layout before returning it.
    for(size_t a=0;a<points.size();++a)for(size_t b=a+1;b<points.size();++b)require((points[a]-points[b]).length()>1e-6,"Optimised graph has coincident vertices; adjust the dimension ranges or graph.");
    zFnGraph(result.graph).create(points,edges);return result;
}
double stairsClearance(const zPointArray& points,const zIntArray& counts,const zIntArray& connects,const zStairsLayoutData& layout,double halfThickness){
    double gap=std::numeric_limits<double>::infinity();intersections(points,counts,connects,layout,[&](auto ids,const auto& samples,const auto& footprint,int owner){double low=1e100,high=-1e100;for(auto weights:samples){double z=0;for(int k=0;k<3;++k)z+=weights[k]*points[ids[k]].z;low=std::min(low,z);high=std::max(high,z);}
        double below=footprint.underside-high-halfThickness;
        double separation=ownedFootprint(layout,footprint,owner)?below:std::max(below,low-halfThickness-footprint.top);
        gap=std::min(gap,separation);
    });require(std::isfinite(gap),"Structural mesh no longer overlaps the walking layout.");return gap;
}
void lowerStairsBelowEnvelope(zPointArray& points,const zIntArray& counts,const zIntArray& connects,const zIntArray& supports,const zStairsLayoutData& layout,double halfThickness,double clearance){
    std::set<int> fixed(supports.begin(),supports.end());
    for(int pass=0;pass<8;++pass){std::vector<double> drops(points.size(),0);double maximum=0;
        intersections(points,counts,connects,layout,[&](auto ids,const auto& samples,const auto& footprint,int owner){if(!ownedFootprint(layout,footprint,owner))return;for(auto weights:samples){double z=0,freeWeight=0;for(int k=0;k<3;++k){z+=weights[k]*points[ids[k]].z;if(!fixed.count(ids[k]))freeWeight+=weights[k];}double violation=z+halfThickness+clearance-footprint.underside;
            if(violation<=1e-6)continue;require(freeWeight>1e-8,"A fixed support intersects a step or landing; edit the source layout.");double drop=(violation+1e-5)/freeWeight;for(int k=0;k<3;++k)if(!fixed.count(ids[k]))drops[ids[k]]=std::max(drops[ids[k]],drop);maximum=std::max(maximum,drop);
        }});
        if(maximum==0)break;for(size_t v=0;v<points.size();++v)points[v].z-=float(drops[v]);
    }require(stairsClearance(points,counts,connects,layout,halfThickness)>=clearance-1e-5,"Structural mesh intersects a different stair flight or cannot clear its step envelope.");
}
}
