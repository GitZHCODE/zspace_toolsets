#include "zTsStairsLayout.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <stdexcept>

namespace zSpace {
namespace {
using V=Eigen::Vector2d;
struct P { V xy; double station; int a,b; double t; };
using Key=std::tuple<int,int,long long,float>;
Key key(const P& p,double z){return {p.a,p.b,llround(p.t*1e8),float(z)};}
P mix(const P& p,const P& q,double t,double station){
    if(t<1e-8){auto r=p;r.station=station;return r;}
    if(t>1-1e-8){auto r=q;r.station=station;return r;}
    std::set<int> ids={p.a,p.b,q.a,q.b};
    if(ids.size()>2)throw std::invalid_argument("Stair face clipping lost its source edge.");
    int a=*ids.begin(),b=*ids.rbegin();
    auto weight=[&](const P& r){return r.a==r.b?double(r.a==b):r.t;};
    return {p.xy+(q.xy-p.xy)*t,station,a,b,weight(p)+(weight(q)-weight(p))*t};
}
std::vector<P> clip(const std::vector<P>& poly,double cut,bool above){
    std::vector<P> out;if(poly.empty())return out;
    for(size_t i=0;i<poly.size();++i){auto p=poly[i],q=poly[(i+1)%poly.size()];
        bool ip=above?p.station>=cut-1e-9:p.station<=cut+1e-9;
        bool iq=above?q.station>=cut-1e-9:q.station<=cut+1e-9;
        if(ip)out.push_back(p);
        if(ip!=iq){double t=(cut-p.station)/(q.station-p.station);out.push_back(mix(p,q,std::clamp(t,0.0,1.0),cut));}
    }return out;
}
double cross(V a,V b){return a.x()*b.y()-a.y()*b.x();}
}
// The walking surface inherits the centre mesh plan and face boundaries. Only
// step contours introduce new cuts; there is no slab extrusion or side closure.
void buildStairsSurface(zStairsLayoutData& layout,const zPointArray& source,
    const zIntArray& counts,const zIntArray& connects,const zParametricStairsSettings& s,bool geometry){
    if(layout.meshFaceFlights.size()!=counts.size())throw std::invalid_argument("Stair face provenance is incomplete.");
    const size_t flights=layout.flightPaths.size();
    std::vector<std::set<int>> owners(source.size());size_t offset=0;
    for(size_t f=0;f<counts.size();++f){for(int i=0;i<counts[f];++i)owners[connects[offset+i]].insert(layout.meshFaceFlights[f]);offset+=counts[f];}
    std::vector<std::map<int,double>> stations(flights);
    auto station=[&](int flight,int v){auto& cache=stations[flight];auto it=cache.find(v);if(it!=cache.end())return it->second;
        const auto& path=layout.flightPaths[flight];double along=0,best=1e100,result=0,total=0;
        V p(source[v].x,source[v].y);
        for(size_t i=1;i<path.size();++i){auto a=layout.graphPoints[path[i-1]],b=layout.graphPoints[path[i]];V x(a.x,a.y),d(b.x-a.x,b.y-a.y);double len=d.norm();double t=std::clamp((p-x).dot(d)/d.squaredNorm(),0.0,1.0);double dist=(p-x-t*d).squaredNorm();if(dist<best){best=dist;result=along+t*len;}along+=len;}total=along;
        // All incident flight patches meet at a common, flat landing seam.
        if(owners[v].size()>1){auto a=layout.graphPoints[path.front()],b=layout.graphPoints[path.back()];result=(p-V(a.x,a.y)).squaredNorm()<=(p-V(b.x,b.y)).squaredNorm()?0:total;}
        cache[v]=result;return result;
    };
    layout.footprints.clear();layout.faceTypes.clear();zPointArray vertices;zIntArray faceCounts,faceConnects;
    std::map<Key,int> welded;std::set<std::pair<Key,Key>> risers;
    auto vertex=[&](const P& p,double z){Key k=key(p,z);auto found=welded.find(k);if(found!=welded.end())return found->second;int id=int(vertices.size());vertices.emplace_back(float(p.xy.x()),float(p.xy.y()),float(z));welded[k]=id;return id;};
    auto horizontal=[&](std::vector<P> poly,double z,int type,int flight){
        if(poly.size()<3)return;std::vector<P> clean;
        for(auto p:poly)if(clean.empty()||key(clean.back(),z)!=key(p,z))clean.push_back(p);
        if(clean.size()>1&&key(clean.front(),z)==key(clean.back(),z))clean.pop_back();
        if(clean.size()<3)return;double area=0;for(size_t i=0;i<clean.size();++i)area+=cross(clean[i].xy,clean[(i+1)%clean.size()].xy);if(std::abs(area)<1e-12)return;if(area<0)std::reverse(clean.begin(),clean.end());
        zPointArray footprint;for(auto p:clean)footprint.emplace_back(float(p.xy.x()),float(p.xy.y()),float(z));
        // Landing faces can be shared by neighbouring flights at the same floor.
        layout.footprints.push_back({footprint,z,z,type==2?-1:flight});
        if(geometry){faceCounts.push_back(int(clean.size()));for(auto p:clean)faceConnects.push_back(vertex(p,z));layout.faceTypes.push_back(type);}
    };
    offset=0;
    for(size_t f=0;f<counts.size();++f){int flight=layout.meshFaceFlights[f];std::vector<P> original;
        for(int i=0;i<counts[f];++i){int v=connects[offset+i];original.push_back({V(source[v].x,source[v].y),station(flight,v),v,v,0});}offset+=counts[f];
        int n=layout.riserCounts[flight];double rise=layout.actualRisers[flight],low=layout.flightFloors[flight][0];
        if(!n){horizontal(original,low,2,flight);continue;}
        double tread=s.flightTreadDepths.empty()?s.treadDepth:s.flightTreadDepths[flight];double first=layout.startReserves[flight]+(layout.usableRuns[flight]-(n-1)*tread)*.5;
        // A saddle scalar field needs a local diagonal so each cut is a single segment.
        bool saddle=false;for(int k=0;k<n;++k){double cut=first+k*tread;int transitions=0;for(size_t i=0;i<original.size();++i)if((original[i].station<cut)!=(original[(i+1)%original.size()].station<cut))++transitions;if(transitions>2)saddle=true;}
        std::vector<std::vector<P>> patches;if(saddle){for(size_t i=1;i+1<original.size();++i)patches.push_back({original[0],original[i],original[i+1]});}else patches.push_back(original);
        for(auto& poly:patches){double minS=1e100,maxS=-1e100;for(auto p:poly){minS=std::min(minS,p.station);maxS=std::max(maxS,p.station);}
            int from=std::clamp(int(std::floor((minS-first)/tread))+1,0,n),to=std::clamp(int(std::floor((maxS-first)/tread))+1,0,n);
            for(int band=from;band<=to;++band){auto part=poly;if(band>0)part=clip(part,first+(band-1)*tread,true);if(band<n)part=clip(part,first+band*tread,false);horizontal(part,low+band*rise,band==0||band==n?2:0,flight);}
            for(int k=0;k<n;++k){double cut=first+k*tread;if(cut<minS-1e-8||cut>maxS+1e-8)continue;std::vector<P> hits;
                auto hit=[&](P p){for(auto q:hits)if(key(p,0)==key(q,0))return;hits.push_back(p);};
                for(size_t i=0;i<poly.size();++i){auto p=poly[i],q=poly[(i+1)%poly.size()];if(std::abs(p.station-cut)<1e-8)hit(p);if((p.station<cut-1e-9&&q.station>cut+1e-9)||(q.station<cut-1e-9&&p.station>cut+1e-9))hit(mix(p,q,(cut-p.station)/(q.station-p.station),cut));}
                if(hits.size()!=2||(hits[0].xy-hits[1].xy).norm()<1e-8)continue;P p=hits[0],q=hits[1];double z0=low+k*rise,z1=z0+rise;
                for(auto r:poly)if(std::abs(r.station-cut)>1e-8){double side=cross(q.xy-p.xy,r.xy-p.xy);if((r.station<cut&&side<0)||(r.station>cut&&side>0))std::swap(p,q);break;}
                auto a=key(p,z0),b=key(q,z0);if(!risers.emplace(std::min(a,b),std::max(a,b)).second)continue;
                V d=q.xy-p.xy,normal(-d.y(),d.x());normal.normalize();normal*=1e-6;zPointArray footprint;for(V v:std::array<V,4>{V(p.xy-normal),V(q.xy-normal),V(q.xy+normal),V(p.xy+normal)})footprint.emplace_back(float(v.x()),float(v.y()),float(z1));layout.footprints.push_back({footprint,z0,z1,flight});
                if(geometry){faceCounts.push_back(4);faceConnects.insert(faceConnects.end(),{vertex(p,z0),vertex(p,z1),vertex(q,z1),vertex(q,z0)});layout.faceTypes.push_back(1);}
            }
        }
    }
    if(geometry){if(faceCounts.empty())throw std::invalid_argument("Stair surface is empty.");zFnMesh(layout.mesh).create(vertices,faceCounts,faceConnects);}
}
}
