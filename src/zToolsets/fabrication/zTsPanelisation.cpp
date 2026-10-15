#include <zspace/zToolsets/fabrication/zTsPanelisation.h>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <climits>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>

namespace zSpace { namespace {
using V = Eigen::Vector3d;
using Edge = std::pair<int,int>;
Edge edge(int a,int b) { return std::minmax(a,b); }
struct DSU {
    std::vector<int> p, rank;
    explicit DSU(int n):p(n),rank(n,0) { std::iota(p.begin(),p.end(),0); }
    int root(int i) { return p[i]==i?i:p[i]=root(p[i]); }
    void join(int a,int b) { a=root(a);b=root(b);if(a==b)return;if(rank[a]<rank[b])std::swap(a,b);p[b]=a;if(rank[a]==rank[b])++rank[a]; }
};
struct Polygon { std::vector<int> vertices; int source; };
struct Triangle { std::array<int,3> vertices; int source, polygon; };
struct Node { V p; int a,b; double t; };
using UV = Eigen::Vector2d;
using Grid = std::pair<int,int>;
struct PanelPlane { V point, normal; };
struct Chart { std::vector<int> faces; std::map<int,UV> uv, flowUV; std::vector<std::pair<double,double>> flowCuts; std::set<double> redCuts; std::map<int,double> bendDepth; int layout=0, source=-1; std::vector<double> stations, crossStations, loopCuts;
    std::vector<V> centreRail, topStarts; std::vector<int> centreCounts; std::vector<std::vector<double>> loopStations; std::vector<std::vector<PanelPlane>> loopPlanes;
    bool flow = false; };
}
void zTsPanelisation::compute(std::vector<zPanelisationSource>& sources,const zPanelisationSettings& s) {
    if(!std::isfinite(s.width)||s.width<=0||!std::isfinite(s.spacing)||s.spacing<=0||!std::isfinite(s.offset)||s.offset<0||s.offset>=s.spacing||
       !std::isfinite(s.seedBand)||s.seedBand<0||!std::isfinite(s.weldTolerance)||s.weldTolerance<0||s.weldTolerance>0.01||
       !std::isfinite(s.creaseAngleDegrees)||s.creaseAngleDegrees<=0||s.creaseAngleDegrees>180||
       !std::isfinite(s.smallPanelArea)||s.smallPanelArea<0||s.axis<0||s.axis>3||sources.empty())
        throw std::invalid_argument("Invalid panel settings: positive spacing and width, 0 <= offset < spacing, axis 0..3, weld tolerance 0..0.01 m.");
    zPanelisationGuides guides;
    std::vector<V> positions; std::vector<Polygon> polygons;
    std::set<Edge> tagged, dividerTags;
    for(int source=0;source<int(sources.size());++source) {
        zPointArray p; zIntArray c,f;zFnMesh fn(sources[source].mesh);fn.getVertexPositions(p);fn.getPolygonData(f,c);
        const auto& options=sources[source];
        if(options.layoutMode<0||options.layoutMode>2||!std::isfinite(options.staggerFraction)||options.staggerFraction<0||options.staggerFraction>=1)throw std::invalid_argument("Invalid panel layout mode or stagger.");
        if(options.layoutMode==1&&(options.alignmentSource<0||options.alignmentSource>=int(sources.size())||options.alignmentDirection.size()!=3))throw std::invalid_argument("Ceiling layout requires an alignment source and direction.");
        if(options.layoutMode==1){V axis(options.alignmentDirection[0],options.alignmentDirection[1],options.alignmentDirection[2]);if(!axis.allFinite()||std::abs(axis.norm()-1)>1e-5)throw std::invalid_argument("Alignment direction must be unit length.");}
        const auto& frame=sources[source].panelFrame;
        if(!frame.empty()) {
            if(frame.size()!=8||!std::all_of(frame.begin(),frame.end(),[](double v){return std::isfinite(v);}))throw std::invalid_argument("Invalid panel layout frame.");
            V u(frame[0],frame[1],frame[2]),v(frame[4],frame[5],frame[6]);
            if(std::abs(u.norm()-1)>1e-6||std::abs(v.norm()-1)>1e-6||std::abs(u.dot(v))>1e-6)throw std::invalid_argument("Panel layout frame axes must be orthonormal.");
        }
        int base=int(positions.size());
        for(auto v:p) { V point(v.x,v.y,v.z);if(!point.allFinite())throw std::invalid_argument("Non-finite input vertex.");positions.push_back(point); }
        size_t offset=0;size_t sourcePolygonStart=polygons.size();
        for(int count:c) {
            if(count<3||count>4||offset+count>f.size())throw std::invalid_argument("Panelisation requires triangles or quads.");
            Polygon polygon;polygon.source=source;
            for(int i=0;i<count;++i) { int v=f[offset++];if(v<0||v>=int(p.size()))throw std::invalid_argument("Invalid mesh index.");polygon.vertices.push_back(base+v); }
            polygons.push_back(std::move(polygon));
        }
        if(offset!=f.size())throw std::invalid_argument("Incomplete polygon connectivity.");
        const auto& ce=sources[source].creaseEdges;
        if(ce.size()%2)throw std::invalid_argument("Invalid crease endpoint pairs.");
        for(size_t i=0;i<ce.size();i+=2) { if(ce[i]<0||ce[i+1]<0||ce[i]>=int(p.size())||ce[i+1]>=int(p.size()))throw std::invalid_argument("Invalid crease vertex.");tagged.insert(edge(base+ce[i],base+ce[i+1])); }
        // A surveyed seed selects a complete existing topology rail. Keep it separate
        // from creases so it partitions panels without restarting their X/Y charts.
        const auto& seeds=sources[source].panelSeamSeeds;
        if(seeds.size()%2)throw std::invalid_argument("Invalid panel seam seed pairs.");
        if(!seeds.empty()) {
            std::map<int,std::set<int>> neighbours;
            std::map<std::pair<int,int>,std::set<int>> turns;
            std::set<int> irregular;
            for(size_t fi=sourcePolygonStart;fi<polygons.size();++fi){const auto& f=polygons[fi].vertices;
                for(size_t j=0;j<f.size();++j){int v=f[j],a=f[(j+f.size()-1)%f.size()],b=f[(j+1)%f.size()];neighbours[v].insert(a);neighbours[v].insert(b);turns[{v,a}].insert(b);turns[{v,b}].insert(a);if(f.size()!=4)irregular.insert(v);}
            }
            auto trace=[&](int a,int b){std::set<Edge> visited;
                while(visited.insert(edge(a,b)).second){dividerTags.insert(edge(a,b));if(neighbours[b].size()!=4||irregular.count(b))break;
                    int next=-1,count=0;for(int v:neighbours[b])if(v!=a&&!turns[{b,a}].count(v)){next=v;++count;}
                    if(count!=1)break;a=b;b=next;
                }
            };
            for(size_t j=0;j<seeds.size();j+=2){int a=seeds[j],b=seeds[j+1];if(a<0||b<0||a>=int(p.size())||b>=int(p.size())||!neighbours[base+a].count(base+b))throw std::invalid_argument("Panel seam seed must be an existing mesh edge.");trace(base+a,base+b);trace(base+b,base+a);}
        }
    }
    if(positions.empty()||polygons.empty()||positions.size()>500000||polygons.size()>500000)throw std::invalid_argument("Empty input or more than 500,000 source vertices/faces.");
    // Stitch only complete, geometrically matching boundary edges. Interior points never weld.
    std::map<Edge,int> incidence;
    std::map<Edge,V> side;
    for(const auto& polygon:polygons){V center=V::Zero();for(int v:polygon.vertices)center+=positions[v];center/=double(polygon.vertices.size());for(size_t i=0;i<polygon.vertices.size();++i){auto e=edge(polygon.vertices[i],polygon.vertices[(i+1)%polygon.vertices.size()]);++incidence[e];V axis=positions[e.second]-positions[e.first];V direction=center-(positions[e.first]+positions[e.second])*0.5;if(axis.squaredNorm()>1e-20)direction-=axis*(direction.dot(axis)/axis.squaredNorm());side[e]=direction.norm()>1e-12?V(direction.normalized()):V::Zero();}}
    for(const auto& e:incidence)if(e.second>2)throw std::invalid_argument("Non-manifold input edge.");
    DSU weld(int(positions.size())),originalComponents(int(positions.size()));int welded=0;
    std::vector<std::array<int,4>> matches;
    for(const auto& item:incidence)originalComponents.join(item.first.first,item.first.second);
    // Coincident vertices at collapsed SubD poles do not establish a unique seam correspondence.
    std::set<int> singularComponents;
    for(const auto& item:incidence)if((positions[item.first.first]-positions[item.first.second]).norm()<1e-10)singularComponents.insert(originalComponents.root(item.first.first));
    if(s.weldTolerance>0) {
        using Cell=std::tuple<long long,long long,long long>;
        std::map<Cell,std::vector<Edge>> buckets;
        auto cell=[&](const V& p){return Cell{(long long)std::floor(p.x()/s.weldTolerance),(long long)std::floor(p.y()/s.weldTolerance),(long long)std::floor(p.z()/s.weldTolerance)};};
        std::set<Edge> paired;
        for(const auto& item:incidence)if(item.second==1) {
            auto e=item.first;V mid=(positions[e.first]+positions[e.second])*0.5;auto [x,y,z]=cell(mid);bool matched=false;
            if(singularComponents.count(originalComponents.root(e.first)))continue;
            for(int dx=-1;dx<=1&&!matched;++dx)for(int dy=-1;dy<=1&&!matched;++dy)for(int dz=-1;dz<=1&&!matched;++dz) {
                auto found=buckets.find({x+dx,y+dy,z+dz});if(found==buckets.end())continue;
                for(auto other:found->second)if(!paired.count(other)&&originalComponents.root(e.first)!=originalComponents.root(other.first)&&side[e].dot(side[other])<0.5) {
                    int a=other.first,b=other.second;
                    if((positions[e.first]-positions[a]).norm()>s.weldTolerance||(positions[e.second]-positions[b]).norm()>s.weldTolerance)std::swap(a,b);
                    if((positions[e.first]-positions[a]).norm()<=s.weldTolerance&&(positions[e.second]-positions[b]).norm()<=s.weldTolerance) {
                        weld.join(e.first,a);weld.join(e.second,b);matches.push_back({e.first,a,e.second,b});paired.insert(other);paired.insert(e);++welded;matched=true;break;
                    }
                }
            }
            buckets[cell(mid)].push_back(e);
        }
    }
    // Reject ambiguous joins locally instead of joining three overlapping skins into one edge.
    for(int pass=0;pass<8;++pass) {
        std::map<Edge,int> counts;
        for(const auto& polygon:polygons) {
            std::vector<int> f;for(int v:polygon.vertices){int root=weld.root(v);if(f.empty()||f.back()!=root)f.push_back(root);}
            if(f.size()>1&&f.front()==f.back())f.pop_back();
            for(size_t j=1;j+1<f.size();++j){int a=f[0],b=f[j],c=f[j+1];if(a==b||b==c||a==c)continue;if((positions[b]-positions[a]).cross(positions[c]-positions[a]).norm()<1e-14)continue;++counts[edge(a,b)];++counts[edge(b,c)];++counts[edge(c,a)];}
        }
        std::set<int> bad;for(const auto& item:counts)if(item.second>2){bad.insert(item.first.first);bad.insert(item.first.second);}
        if(bad.empty())break;
        std::vector<std::array<int,4>> valid;for(auto match:matches){bool keep=true;for(int v:match)if(bad.count(weld.root(v)))keep=false;if(keep)valid.push_back(match);}
        if(valid.size()==matches.size())throw std::invalid_argument("Input triangulation is non-manifold.");
        matches=std::move(valid);weld=DSU(int(positions.size()));for(auto match:matches){weld.join(match[0],match[1]);weld.join(match[2],match[3]);}welded=int(matches.size());
    }
    std::vector<V> p;std::map<int,int> remap;std::vector<int> ids(positions.size());
    for(int i=0;i<int(positions.size());++i) { int root=weld.root(i);auto it=remap.emplace(root,int(p.size()));if(it.second)p.push_back(positions[root]);ids[i]=it.first->second; }
    std::set<Edge> creases;for(auto e:tagged)creases.insert(edge(ids[e.first],ids[e.second]));
    std::set<Edge> dividers;for(auto e:dividerTags)dividers.insert(edge(ids[e.first],ids[e.second]));
    std::vector<Triangle> triangles;std::map<Edge,std::vector<int>> edges;
    for(int polygonId=0;polygonId<int(polygons.size());++polygonId) {
        auto& polygon=polygons[polygonId];
        for(int& v:polygon.vertices)v=ids[v];
        polygon.vertices.erase(std::unique(polygon.vertices.begin(),polygon.vertices.end()),polygon.vertices.end());
        if(polygon.vertices.size()>1&&polygon.vertices.front()==polygon.vertices.back())polygon.vertices.pop_back();
        if(polygon.vertices.size()<3)continue;
        if(std::set<int>(polygon.vertices.begin(),polygon.vertices.end()).size()!=polygon.vertices.size())throw std::invalid_argument("Stitching collapsed a non-adjacent face vertex; reduce weld tolerance.");
        for(size_t j=1;j+1<polygon.vertices.size();++j) {
            Triangle t{{polygon.vertices[0],polygon.vertices[j],polygon.vertices[j+1]},polygon.source,polygonId};
            // Saved Rhino meshes can contain zero-area triangles at collapsed surface poles.
            if((p[t.vertices[1]]-p[t.vertices[0]]).cross(p[t.vertices[2]]-p[t.vertices[0]]).norm()<1e-14)continue;
            int index=int(triangles.size());triangles.push_back(t);
            for(int k=0;k<3;++k)edges[edge(t.vertices[k],t.vertices[(k+1)%3])].push_back(index);
        }
    }
    std::vector<bool> activePolygon(polygons.size(),false);for(const auto& t:triangles)activePolygon[t.polygon]=true;
    for(size_t i=0;i<polygons.size();++i)if(!activePolygon[i])polygons[i].vertices.clear();
    auto normal=[&](const Triangle& t){return (p[t.vertices[1]]-p[t.vertices[0]]).cross(p[t.vertices[2]]-p[t.vertices[0]]).normalized();};
    std::vector<bool> boundary(p.size(),false);DSU components(int(p.size()));
    std::vector<std::vector<std::pair<int,double>>> adjacency(p.size());
    for(const auto& item:edges) {
        auto [a,b]=item.first;const auto& faces=item.second;
        if(faces.size()>2)throw std::invalid_argument("Stitching produced a non-manifold edge at source "+std::to_string(triangles[faces[0]].source)+"; reduce weld tolerance.");
        if(faces.size()==1)boundary[a]=boundary[b]=true;
        if(faces.size()==2&&s.useCreases) {
            double cosine=std::clamp(normal(triangles[faces[0]]).dot(normal(triangles[faces[1]])),-1.0,1.0);
            if(std::acos(cosine)*180.0/3.141592653589793>=s.creaseAngleDegrees)creases.insert(item.first);
        }
        double length=(p[a]-p[b]).norm();components.join(a,b);adjacency[a].push_back({b,length});adjacency[b].push_back({a,length});
    }
    if(!s.useCreases)creases.clear();
    std::map<int,std::vector<int>> groups;for(int i=0;i<int(p.size());++i)if(!adjacency[i].empty())groups[components.root(i)].push_back(i);
    // Grow intrinsic lattice charts through opposite quad edges. The grid rotates with
    // the surface, so walls, soffits and curved canopies use the same construction.
    std::map<Edge,std::vector<int>> polygonEdges;
    for(int i=0;i<int(polygons.size());++i){const auto& f=polygons[i].vertices;for(size_t j=0;j<f.size();++j)polygonEdges[edge(f[j],f[(j+1)%f.size()])].push_back(i);}
    DSU smooth(int(polygons.size()));
    for(const auto& item:polygonEdges)if(item.second.size()==2&&!creases.count(item.first))smooth.join(item.second[0],item.second[1]);
    std::set<int> mixed;for(int i=0;i<int(polygons.size());++i)if(polygons[i].vertices.size()==3)mixed.insert(smooth.root(i));
    std::vector<int> chartIds(polygons.size(),-1);std::vector<Chart> charts;
    for(int seed=0;seed<int(polygons.size());++seed)if(chartIds[seed]<0&&polygons[seed].vertices.size()>=3) {
        int chartId=int(charts.size());Chart chart;chart.source=polygons[seed].source;chart.layout=sources[chart.source].layoutMode;
        if(polygons[seed].vertices.size()==4&&!mixed.count(smooth.root(seed))) {
            std::map<int,Grid> grid;std::map<Grid,int> occupied;
            const Grid initial[]={{0,0},{1,0},{1,1},{0,1}};
            for(int j=0;j<4;++j){int v=polygons[seed].vertices[j];grid[v]=initial[j];occupied[initial[j]]=v;}
            chartIds[seed]=chartId;chart.faces.push_back(seed);
            for(size_t cursor=0;cursor<chart.faces.size();++cursor) {
                int parent=chart.faces[cursor];const auto& f=polygons[parent].vertices;
                for(int j=0;j<4;++j){int a=f[j],b=f[(j+1)%4];auto e=edge(a,b);if(creases.count(e)&&chart.layout==0)continue;
                    for(int next:polygonEdges[e])if(chartIds[next]<0&&polygons[next].vertices.size()==4&&((chart.layout==0&&sources[polygons[next].source].layoutMode==0)||(chart.layout!=0&&polygons[next].source==chart.source))){
                        const auto& n=polygons[next].vertices;int k=0;while(k<4&&n[k]!=a)++k;
                        if(k==4)continue;int step=n[(k+1)%4]==b?1:n[(k+3)%4]==b?-1:0;if(!step)continue;
                        auto ga=grid.at(a),gb=grid.at(b),inside=grid.at(f[(j+3)%4]);
                        Grid outward={ga.first-inside.first,ga.second-inside.second};
                        std::array<std::pair<int,Grid>,2> proposed={{{n[(k+2*step+8)%4],{gb.first+outward.first,gb.second+outward.second}},
                            {n[(k+3*step+12)%4],{ga.first+outward.first,ga.second+outward.second}}}};
                        bool valid=true;
                        for(const auto& item:proposed){auto known=grid.find(item.first);auto cell=occupied.find(item.second);if((known!=grid.end()&&known->second!=item.second)||(cell!=occupied.end()&&cell->second!=item.first))valid=false;}
                        if(!valid)continue; // extraordinary/periodic topology starts another chart
                        for(const auto& item:proposed){grid[item.first]=item.second;occupied[item.second]=item.first;}
                        chartIds[next]=chartId;chart.faces.push_back(next);
                    }
                }
            }
            // Review-only flow guides: equal arc-length seeds on the highest chart side,
            // transported through opposite quad edges. Refinement is sampled at every row.
            int guideSource=polygons[seed].source;
            if(sources[guideSource].topologyGuides) {
                int minimum[2]={INT_MAX,INT_MAX},maximum[2]={INT_MIN,INT_MIN};
                for(auto item:grid){minimum[0]=std::min(minimum[0],item.second.first);maximum[0]=std::max(maximum[0],item.second.first);minimum[1]=std::min(minimum[1],item.second.second);maximum[1]=std::max(maximum[1],item.second.second);}
                int down=0,top=minimum[0],sign=1;double best=-std::numeric_limits<double>::infinity();
                for(int d=0;d<2;++d)for(int side:{minimum[d],maximum[d]}){double height=0;int count=0;for(auto item:grid)if((d?item.second.second:item.second.first)==side){height+=p[item.first].z();++count;}if(count>1&&height/count>best){best=height/count;down=d;top=side;sign=side==minimum[d]?1:-1;}}
                using Row=std::map<int,V>;std::map<int,Row> rows;
                for(auto item:grid){int depth=down?item.second.second:item.second.first,across=down?item.second.first:item.second.second;rows[(depth-top)*sign][across]=p[item.first];}
                std::map<int,std::map<int,double>> pathDepths;std::map<int,double> depths;
                const Row* lastRow=nullptr;int lastDepth=0;
                for(const auto& row:rows){double sum=0;for(const auto& q:row.second){double length=0;if(lastRow){auto previous=lastRow->find(q.first);if(previous!=lastRow->end())length=pathDepths.at(lastDepth).at(q.first)+(q.second-previous->second).norm();}pathDepths[row.first][q.first]=length;sum+=length;}depths[row.first]=sum/row.second.size();lastRow=&row.second;lastDepth=row.first;}
                // The straight wall rows end where a column turns away from its
                // initial tangent. Ten degrees ignores small limit-mesh noise.
                if(!sources[guideSource].panelSeamSeeds.empty())for(const auto& topPoint:rows.begin()->second){
                    int u=topPoint.first;V previous=topPoint.second,initial=V::Zero();double previousLength=0;
                    chart.bendDepth[u]=1e30;
                    for(auto ri=std::next(rows.begin());ri!=rows.end();++ri){auto q=ri->second.find(u);if(q==ri->second.end())break;V tangent=q->second-previous;
                        if(tangent.norm()>1e-10){tangent.normalize();if(initial.norm()<0.5)initial=tangent;
                            else if(initial.dot(tangent)<std::cos(10.0*3.141592653589793/180.0)){chart.bendDepth[u]=previousLength;break;}}
                        previous=q->second;previousLength=pathDepths.at(ri->first).at(u);
                    }
                }
                auto pathDepth=[&](int row,double u){const auto& samples=pathDepths.at(row);auto b=samples.lower_bound(int(std::ceil(u)));if(b==samples.end())return samples.rbegin()->second;if(b==samples.begin()||std::abs(b->first-u)<1e-10)return b->second;auto a=std::prev(b);return a->second+(u-a->first)/(b->first-a->first)*(b->second-a->second);};
                for(auto item:grid){int depth=((down?item.second.second:item.second.first)-top)*sign,across=down?item.second.first:item.second.second;chart.flowUV[item.first]=UV(across,pathDepths.at(depth).at(across));}
                auto point=[](const Row& row,double u,V& out){auto b=row.lower_bound(int(std::ceil(u)));if(b==row.end())return false;if(std::abs(b->first-u)<1e-12){out=b->second;return true;}if(b==row.begin())return false;auto a=std::prev(b);if(b->first-a->first!=1||u<a->first||u>b->first)return false;out=a->second+(u-a->first)*(b->second-a->second);return true;};
                auto width=[&](const Row& row,double a,double b){V previous;if(!point(row,a,previous))return -1.0;double length=0;for(int k=int(std::floor(a))+1;k<b-1e-9;++k){V next;if(!point(row,k,next))return -1.0;length+=(next-previous).norm();previous=next;}V last;if(!point(row,b,last))return -1.0;return length+(last-previous).norm();};
                auto station=[&](const Row& row,double a,double b,double distance){double low=a,high=b;for(int k=0;k<36;++k){double mid=(low+high)/2,w=width(row,a,mid);if(w<0)return a;if(w<distance)low=mid;else high=mid;}double u=(low+high)/2;return std::abs(u-std::round(u))<1e-7?std::round(u):u;};
                auto segment=[&](const V& a,const V& b,int kind){if((b-a).norm()<1e-10)return;for(auto q:{a,b})guides.positions.insert(guides.positions.end(),{q.x(),q.y(),q.z()});guides.kinds.push_back(kind);guides.sources.push_back(guideSource);};
                struct Path {double u;int kind;V last;};std::vector<Path> paths;const auto& first=rows.begin()->second;double a=first.begin()->first,b=first.rbegin()->first,total=width(first,a,b);
                if(total<=0){++guides.unsupportedCharts;}else {
                    // Creases entering the wall are mandatory top anchors. A crease
                    // running along the top edge is not a station at every mesh vertex.
                    // Redistribute each intervening span, avoiding a second equal-grid
                    // cut just beside an authored crease (including partial creases).
                    std::set<double> anchors={a,b};
                    for(auto e:creases){auto ia=grid.find(e.first),ib=grid.find(e.second);if(ia==grid.end()||ib==grid.end())continue;
                        int da=down?ia->second.second:ia->second.first,db=down?ib->second.second:ib->second.first;
                        int ua=down?ia->second.first:ia->second.second,ub=down?ib->second.first:ib->second.second;
                        if(ua==ub&&da!=db&&(da==top||db==top)&&ua>=a&&ua<=b){anchors.insert(ua);chart.redCuts.insert(ua);}
                    }
                    std::vector<double> starts={a};
                    for(auto it=std::next(anchors.begin());it!=anchors.end();++it){double left=*std::prev(it),right=*it,length=width(first,left,right);
                        int intervals=std::max(1,int(std::ceil(length/s.spacing-1e-6)));if(starts.size()+intervals>10001)throw std::invalid_argument("Too many topology guide seeds.");
                        for(int i=1;i<=intervals;++i)starts.push_back(i==intervals?right:station(first,left,right,length*i/intervals));
                    }
                    for(double u:starts){V q;if(point(first,u,q)){paths.push_back({u,0,q});chart.topStarts.push_back(q);chart.flowCuts.push_back({u,0});++guides.seeds;}}
                    chart.flow=true;
                    const Row* previous=&first;int previousRow=rows.begin()->first;
                    for(auto rowIt=std::next(rows.begin());rowIt!=rows.end();++rowIt){const auto& row=rowIt->second;std::vector<Path> active;
                        for(auto path:paths){V q;if(point(row,path.u,q)){segment(path.last,q,path.kind);path.last=q;active.push_back(path);}}
                        std::vector<Path> added;
                        for(size_t i=1;i<active.size();++i){double left=active[i-1].u,right=active[i].u,w=width(row,left,right);if(w<=0)continue;
                            int divisions=w>s.spacing*(1+1e-6)?int(std::ceil(w/s.spacing)):1;
                            if(divisions>10000)throw std::invalid_argument("Too many adaptive topology guides.");
                            guides.maxSampledWidth=std::max(guides.maxSampledWidth,w/divisions);
                            for(int j=1;j<divisions;++j){double u=station(row,left,right,w*j/divisions);V q,old;point(row,u,q);double birth=pathDepth(rowIt->first,u);
                                if(point(*previous,u,old)){double oldWidth=width(*previous,left,right),fraction=oldWidth>=0&&w>oldWidth?std::clamp((s.spacing-oldWidth)/(w-oldWidth),0.0,1.0):1.0;segment(old+fraction*(q-old),q,1);birth=pathDepth(previousRow,u)+fraction*(pathDepth(rowIt->first,u)-pathDepth(previousRow,u));}
                                added.push_back({u,1,q});chart.flowCuts.push_back({u,birth});++guides.branches;
                            }
                        }
                        active.insert(active.end(),added.begin(),added.end());std::sort(active.begin(),active.end(),[](const Path& a,const Path& b){return a.u<b.u;});paths=std::move(active);previous=&row;previousRow=rowIt->first;
                        if(guides.kinds.size()>1000000)throw std::invalid_argument("Too many topology guide segments.");
                    }
                }
            }
            // Integrate physical edge lengths per lattice interval, not world projections.
            std::map<int,std::pair<double,int>> lengths[2];double creaseScore[2]={0,0},extent[2]={0,0};std::set<Edge> visited;
            for(int f:chart.faces){const auto& poly=polygons[f].vertices;for(int j=0;j<4;++j){int a=poly[j],b=poly[(j+1)%4];auto e=edge(a,b);if(!visited.insert(e).second)continue;
                auto ga=grid.at(a),gb=grid.at(b);int direction=ga.first!=gb.first?0:1;int index=direction==0?std::min(ga.first,gb.first):std::min(ga.second,gb.second);
                double length=(p[a]-p[b]).norm();auto& interval=lengths[direction][index];interval.first+=length;++interval.second;
                if(creases.count(e))creaseScore[direction]+=length;
            }}
            std::map<int,double> coordinates[2];
            for(int d=0;d<2;++d){double sum=0;for(auto item:lengths[d]){coordinates[d][item.first]=sum;sum+=item.second.first/item.second.second;coordinates[d][item.first+1]=sum;}extent[d]=sum;}
            // X runs along crease rails; on an uncreased patch choose its longer row family.
            int along=std::max(creaseScore[0],creaseScore[1])>1e-8?(creaseScore[1]>creaseScore[0]?1:0):(extent[1]>extent[0]?1:0);
            for(auto item:grid){UV q(coordinates[0].at(item.second.first),coordinates[1].at(item.second.second));if(along==1)std::swap(q.x(),q.y());chart.uv[item.first]=q;}
            if(chart.layout){
                std::map<int,std::pair<V,int>> centres;
                for(auto item:grid){int u=along?item.second.second:item.second.first;auto found=centres.emplace(u,std::make_pair(V::Zero().eval(),0));found.first->second.first+=p[item.first];++found.first->second.second;}
                for(auto item:centres){chart.stations.push_back(coordinates[along].at(item.first));chart.centreRail.push_back(item.second.first/item.second.second);chart.centreCounts.push_back(item.second.second);}
                for(auto item:coordinates[1-along])chart.crossStations.push_back(item.second);
                std::set<double> loops={0,extent[1-along]};
                for(auto e:creases)if(chart.uv.count(e.first)&&chart.uv.count(e.second)){auto a=chart.uv.at(e.first),b=chart.uv.at(e.second);if(std::abs(a.y()-b.y())<1e-8&&std::abs(a.x()-b.x())>1e-8)loops.insert(a.y());}
                chart.loopCuts.assign(loops.begin(),loops.end());
            }
        } else {
            if(sources[polygons[seed].source].topologyGuides)++guides.unsupportedCharts;
            // Triangle-only patches use a conformal intrinsic chart, with two boundary pins.
            chartIds[seed]=chartId;chart.faces.push_back(seed);
            auto polygonNormal=[&](int face){const auto& v=polygons[face].vertices;V n=(p[v[1]]-p[v[0]]).cross(p[v[2]]-p[v[0]]);return n.norm()>1e-12?V(n.normalized()):V::Zero();};
            V seedNormal=polygonNormal(seed);
            for(size_t cursor=0;cursor<chart.faces.size();++cursor){const auto& f=polygons[chart.faces[cursor]].vertices;for(size_t j=0;j<f.size();++j){auto e=edge(f[j],f[(j+1)%f.size()]);if(creases.count(e))continue;for(int next:polygonEdges[e])if(chartIds[next]<0&&polygons[next].vertices.size()>=3&&polygonNormal(next).dot(seedNormal)>0.1&&std::count_if(polygons[next].vertices.begin(),polygons[next].vertices.end(),[&](int v){return std::find(f.begin(),f.end(),v)!=f.end();})==2){chartIds[next]=chartId;chart.faces.push_back(next);}}}
            std::vector<std::array<int,3>> chartTriangles;
            for(int f:chart.faces){const auto& vs=polygons[f].vertices;for(size_t j=1;j+1<vs.size();++j)if((p[vs[j]]-p[vs[0]]).cross(p[vs[j+1]]-p[vs[0]]).norm()>1e-14)chartTriangles.push_back({vs[0],vs[j],vs[j+1]});}
            std::map<int,int> local;std::map<Edge,int> ce;
            for(int f:chart.faces){const auto& vs=polygons[f].vertices;for(int j=0;j<int(vs.size());++j)++ce[edge(vs[j],vs[(j+1)%vs.size()])];}
            for(const auto& t:chartTriangles)for(int v:t)local.emplace(v,int(local.size()));
            std::vector<int> boundaryIds;for(auto item:ce)if(item.second==1){if(local.count(item.first.first))boundaryIds.push_back(item.first.first);if(local.count(item.first.second))boundaryIds.push_back(item.first.second);}
            if(boundaryIds.empty())throw std::invalid_argument("Closed triangle patch at source "+std::to_string(polygons[seed].source)+" seed "+std::to_string(seed)+" faces "+std::to_string(chart.faces.size()));
            int anchor=boundaryIds[0],other=anchor;double span=0;
            for(int pass=0;pass<2;++pass){span=0;for(int v:boundaryIds){double d=(p[v]-p[anchor]).squaredNorm();if(d>span){span=d;other=v;}}if(pass==0)anchor=other;}
            span=std::sqrt(span);if(span<1e-10)throw std::invalid_argument("Degenerate triangle chart.");
            int n=int(local.size());std::map<int,double> fixed={{2*local.at(anchor),0},{2*local.at(anchor)+1,0},{2*local.at(other),span},{2*local.at(other)+1,0}};
            std::vector<int> freeIndex(2*n,-1);int count=0;for(int i=0;i<2*n;++i)if(!fixed.count(i))freeIndex[i]=count++;
            std::vector<Eigen::Triplet<double>> entries;Eigen::VectorXd rhs=Eigen::VectorXd::Zero(2*chartTriangles.size());int row=0;double surfaceArea=0;
            for(const auto& vs:chartTriangles){V e=p[vs[1]]-p[vs[0]],r=p[vs[2]]-p[vs[0]];double len=e.norm(),x=r.dot(e)/len,y=e.cross(r).norm()/len;
                if(y<1e-20||len<1e-20){row+=2;continue;}double area=len*y/2;surfaceArea+=area;double w=1/std::sqrt(2*len*y);
                UV gradients[]={{-y,x-len},{y,-x},{0,len}};
                for(int j=0;j<3;++j){int u=2*local.at(vs[j]),v=u+1;double gx=gradients[j].x()*w,gy=gradients[j].y()*w;
                    for(auto coefficient:std::array<std::tuple<int,int,double>,4>{{{row,u,gx},{row,v,-gy},{row+1,u,gy},{row+1,v,gx}}}){auto [r0,c,value]=coefficient;if(fixed.count(c))rhs[r0]-=value*fixed.at(c);else entries.emplace_back(r0,freeIndex[c],value);}}
                row+=2;
            }
            Eigen::SparseMatrix<double> a(rhs.size(),count);a.setFromTriplets(entries.begin(),entries.end());Eigen::SparseMatrix<double> system=a.transpose()*a;
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;solver.compute(system);if(solver.info()!=Eigen::Success)throw std::invalid_argument("Cannot parameterise triangle patch at source "+std::to_string(polygons[seed].source)+" faces "+std::to_string(chart.faces.size()));
            Eigen::VectorXd solution=solver.solve(a.transpose()*rhs);if(solver.info()!=Eigen::Success||!solution.allFinite())throw std::invalid_argument("Triangle chart solve failed.");
            for(auto item:local){UV q;for(int d=0;d<2;++d){int c=2*item.second+d;q[d]=fixed.count(c)?fixed.at(c):solution[freeIndex[c]];}chart.uv[item.first]=q;}
            double flatArea=0;for(const auto& v:chartTriangles){UV a=chart.uv.at(v[1])-chart.uv.at(v[0]),b=chart.uv.at(v[2])-chart.uv.at(v[0]);flatArea+=std::abs(a.x()*b.y()-a.y()*b.x())/2;}
            if(flatArea<1e-12)throw std::invalid_argument("Collapsed triangle chart; add a seam.");
            // Align to boundary rails rather than the diagonal joining the two solver pins.
            double cx=0,cy=0;for(auto item:ce)if(item.second==1&&chart.uv.count(item.first.first)&&chart.uv.count(item.first.second)){UV d=chart.uv.at(item.first.second)-chart.uv.at(item.first.first);double angle=std::atan2(d.y(),d.x()),w=d.norm();cx+=w*std::cos(4*angle);cy+=w*std::sin(4*angle);}
            double angle=std::atan2(cy,cx)/4,scale=std::sqrt(surfaceArea/flatArea);UV minimum=UV::Constant(std::numeric_limits<double>::infinity());
            for(auto& item:chart.uv){auto q=item.second;item.second=scale*UV(std::cos(angle)*q.x()+std::sin(angle)*q.y(),-std::sin(angle)*q.x()+std::cos(angle)*q.y());minimum=minimum.cwiseMin(item.second);}
            for(auto& item:chart.uv)item.second-=minimum;
        }
        if(chart.flow)chart.uv=chart.flowUV;
        charts.push_back(std::move(chart));
    }
    struct CeilingRail { std::vector<V> points; std::vector<double> distances; };
    auto extractRail=[&](const Chart& chart,double cross){
        CeilingRail rail;std::map<double,V> ordered;
        int full=*std::max_element(chart.centreCounts.begin(),chart.centreCounts.end());
        double first=1e30,last=-1e30;
        for(size_t i=0;i<chart.stations.size();++i)if(chart.centreCounts[i]>=0.9*full){first=std::min(first,chart.stations[i]);last=std::max(last,chart.stations[i]);}
        // Partial corner returns can carry an outer chart coordinate that does
        // not extend along the roof. Select the nearest complete native rail.
        std::map<double,std::map<double,V>> candidates;
        for(auto item:chart.uv)if(item.second.x()>=first-1e-7&&item.second.x()<=last+1e-7)candidates[item.second.y()][item.second.x()]=p[item.first];
        size_t longest=0;for(const auto& candidate:candidates)longest=std::max(longest,candidate.second.size());
        double best=1e30;for(const auto& candidate:candidates)if(candidate.second.size()>=std::max<size_t>(2,size_t(0.9*longest))&&std::abs(candidate.first-cross)<best){best=std::abs(candidate.first-cross);ordered=candidate.second;}

        for(auto item:ordered){double d=rail.points.empty()?0:rail.distances.back()+(item.second-rail.points.back()).norm();rail.points.push_back(item.second);rail.distances.push_back(d);}
        if(rail.points.size()<2||rail.distances.back()<1e-8)throw std::invalid_argument("Cannot extract ceiling boundary rail.");
        return rail;
    };
    auto railPoint=[](const CeilingRail& rail,double fraction)->V{
        double d=std::clamp(fraction,0.0,1.0)*rail.distances.back();auto it=std::lower_bound(rail.distances.begin(),rail.distances.end(),d);
        size_t i=std::clamp<size_t>(it-rail.distances.begin(),1,rail.distances.size()-1);double span=rail.distances[i]-rail.distances[i-1];
        return rail.points[i-1]+(span>1e-12?(d-rail.distances[i-1])/span:0)*(rail.points[i]-rail.points[i-1]);
    };
    int referenceId=-1;for(int i=0;i<int(charts.size());++i)if(charts[i].layout==1&&(referenceId<0||charts[i].faces.size()>charts[referenceId].faces.size()))referenceId=i;
    int ceilingDivisions=1;
    if(referenceId>=0){const auto& reference=charts[referenceId];double shortest=1e30;
        for(double cross:reference.loopCuts)shortest=std::min(shortest,extractRail(reference,cross).distances.back());
        ceilingDivisions=std::max(1,int(std::ceil(shortest/(2*s.spacing))));
        if(ceilingDivisions>10000)throw std::invalid_argument("Too many ceiling divisions.");
    }
    int archReferenceId=-1;for(int i=0;i<int(charts.size());++i)if(charts[i].layout==2&&(archReferenceId<0||charts[i].faces.size()>charts[archReferenceId].faces.size()))archReferenceId=i;
    int archDivisions=1;
    if(archReferenceId>=0){const auto& reference=charts[archReferenceId];double shortest=1e30;
        for(double cross:reference.loopCuts)shortest=std::min(shortest,extractRail(reference,cross).distances.back());
        archDivisions=std::max(1,int(std::ceil(shortest/(2*s.spacing))));if(archDivisions%2==0)++archDivisions;
        if(archDivisions>10000)throw std::invalid_argument("Too many arch divisions.");
    }
    auto nearest=[](const std::vector<double>& values,double target){auto it=std::lower_bound(values.begin(),values.end(),target);if(it==values.end())return values.back();if(it!=values.begin()&&target-*std::prev(it)<=*it-target)return *std::prev(it);return *it;};
    for(auto& chart:charts)if(chart.layout){
        if(chart.stations.size()<2||chart.loopCuts.size()<2)throw std::invalid_argument("Custom panel layout requires a regular quad strip.");
        const auto& options=sources[chart.source];
        if(chart.layout==1){std::set<double> loops(chart.loopCuts.begin(),chart.loopCuts.end());
            for(size_t i=1;i<chart.loopCuts.size();++i){double a=chart.loopCuts[i-1],b=chart.loopCuts[i];int n=std::max(1,int(std::ceil((b-a)/s.width-1e-6)));if(n>10000)throw std::invalid_argument("Too many ceiling loops.");for(int j=1;j<n;++j)loops.insert(nearest(chart.crossStations,a+(b-a)*j/n));}
            chart.loopCuts.assign(loops.begin(),loops.end());
        }
        double widestLoop=0;for(size_t i=1;i<chart.loopCuts.size();++i)widestLoop=std::max(widestLoop,chart.loopCuts[i]-chart.loopCuts[i-1]);
        int panelLoop=0;
        for(size_t loop=0;loop+1<chart.loopCuts.size();++loop){double v0=chart.loopCuts[loop],v1=chart.loopCuts[loop+1];
            {
                const auto left=extractRail(chart,v0),right=extractRail(chart,v1);
                std::vector<PanelPlane> planes;std::vector<double> cuts;
                int parity=panelLoop%2;
                int divisions=chart.layout==1?ceilingDivisions:archDivisions;
                for(int j=parity?1:2;j<divisions;j+=2){double fraction=double(j)/divisions;
                    V a=railPoint(left,fraction),b=railPoint(right,fraction),across=b-a;
                    V leftForward=railPoint(left,std::min(1.0,fraction+0.001))-railPoint(left,std::max(0.0,fraction-0.001));
                    V rightForward=railPoint(right,std::min(1.0,fraction+0.001))-railPoint(right,std::max(0.0,fraction-0.001));
                    V forward=leftForward+rightForward;
                    V normal=forward;if(across.squaredNorm()>1e-16)normal-=across*forward.dot(across)/across.squaredNorm();
                    V origin=a;
                    if(chart.layout==2&&options.alignmentDirection.size()==3){
                        V axis(options.alignmentDirection[0],options.alignmentDirection[1],0);axis.normalize();
                        double horizontal=across.dot(axis);
                        normal=axis*(-across.z())+V(0,0,horizontal);
                        if(normal.norm()<1e-8)normal=axis*forward.dot(axis)+V(0,0,forward.z());
                    }
                    // Equal fractions of unequal rails can skew a straight reach.
                    // Lock its section normal to the surveyed reach axis; retain
                    // the rail-to-rail fan only where the rails actually turn.
                    if(options.alignmentDirection.size()==3&&leftForward.norm()>1e-10&&rightForward.norm()>1e-10){
                        V axis(options.alignmentDirection[0],options.alignmentDirection[1],0);axis.normalize();
                        std::array<V,3> reachAxes={axis,V(-axis.y(),axis.x(),0),V(0,0,1)};
                        for(const V& reach:reachAxes)if(std::abs(leftForward.normalized().dot(reach))>0.984807753&&std::abs(rightForward.normalized().dot(reach))>0.984807753){
                            normal=reach;origin=(a+b)*0.5;break;
                        }
                    }
                    if(normal.dot(forward)<0)normal=-normal;
                    if(normal.norm()<1e-10)continue;normal.normalize();planes.push_back({origin,normal});cuts.push_back(fraction);
                }
                if(v1-v0>=widestLoop*0.2)++panelLoop;
                chart.loopPlanes.push_back(std::move(planes));chart.loopStations.push_back(std::move(cuts));continue;
            }
        }
    }
    // Extraordinary vertices leave small end patches outside the main quad
    // chart. They inherit its stations instead of repeating the full n divisions.
    for(int mainId:std::array<int,2>{referenceId,archReferenceId})if(mainId>=0){const auto& reference=charts[mainId];
        std::vector<std::vector<V>> bandSamples;
        for(size_t band=0;band+1<reference.loopCuts.size();++band){auto a=extractRail(reference,reference.loopCuts[band]),b=extractRail(reference,reference.loopCuts[band+1]);std::vector<V> samples;
            for(int j=0;j<=128;++j)samples.push_back((railPoint(a,j/128.0)+railPoint(b,j/128.0))*0.5);bandSamples.push_back(std::move(samples));}
        for(size_t ci=0;ci<charts.size();++ci){auto& chart=charts[ci];if(int(ci)==mainId||chart.layout!=reference.layout||chart.source!=reference.source)continue;
            for(size_t band=0;band<chart.loopPlanes.size();++band){V centre=V::Zero();int count=0;
                for(auto item:chart.uv)if(item.second.y()>=chart.loopCuts[band]-1e-8&&item.second.y()<=chart.loopCuts[band+1]+1e-8){centre+=p[item.first];++count;}
                if(!count)continue;centre/=count;double best=1e30;size_t nearestBand=0;
                for(size_t k=0;k<bandSamples.size();++k)for(const V& point:bandSamples[k]){double distance=(point-centre).squaredNorm();if(distance<best){best=distance;nearestBand=k;}}
                chart.loopPlanes[band]=reference.loopPlanes[nearestBand];chart.loopStations[band]=reference.loopStations[nearestBand];
            }
        }
    }
    // Clip in both intrinsic coordinates. Edge parameters identify shared intersections,
    // independent of chart direction, source object, or clipping order.
    std::vector<Node> nodes;for(int i=0;i<int(p.size());++i)nodes.push_back({p[i],i,i,0});
    std::map<Edge,std::map<long long,int>> edgeCuts;
    std::map<std::tuple<int,long long,long long>,int> insideCuts;
    struct Clip { UV q; V bary; };
    auto nodeId=[&](const Clip& v,int ti){const auto& t=triangles[ti];for(int j=0;j<3;++j)if(v.bary[j]>1-1e-9)return t.vertices[j];
        if(charts[chartIds[t.polygon]].layout){V position=v.bary[0]*p[t.vertices[0]]+v.bary[1]*p[t.vertices[1]]+v.bary[2]*p[t.vertices[2]];
            for(int j=0;j<3;++j)if((position-p[t.vertices[j]]).squaredNorm()<1e-14)return t.vertices[j];}

        for(int j=0;j<3;++j)if(std::abs(v.bary[j])<1e-9){int i=(j+1)%3,k=(j+2)%3,a=t.vertices[i],b=t.vertices[k];double parameter=v.bary[k]/(v.bary[i]+v.bary[k]);if(a>b){std::swap(a,b);parameter=1-parameter;}
            auto key=(long long)std::llround(parameter*1e9);auto& cache=edgeCuts[{a,b}];auto found=cache.find(key);if(found!=cache.end())return found->second;
            int id=int(nodes.size());nodes.push_back({p[a]+parameter*(p[b]-p[a]),a,b,parameter});cache[key]=id;return id;}
        auto key=std::make_tuple(ti,(long long)std::llround(v.bary[0]*1e9),(long long)std::llround(v.bary[1]*1e9));auto found=insideCuts.find(key);if(found!=insideCuts.end())return found->second;
        int id=int(nodes.size());nodes.push_back({v.bary[0]*p[t.vertices[0]]+v.bary[1]*p[t.vertices[1]]+v.bary[2]*p[t.vertices[2]],-1,-1,0});insideCuts[key]=id;return id;
    };
    auto clip=[](std::vector<Clip> f,int dimension,double level,bool above){std::vector<Clip> out;if(f.empty())return out;
        for(size_t i=0;i<f.size();++i){const auto& a=f[i];const auto& b=f[(i+1)%f.size()];double da=a.q[dimension]-level,db=b.q[dimension]-level;bool ina=above?da>=-1e-10:da<=1e-10,inb=above?db>=-1e-10:db<=1e-10;if(ina)out.push_back(a);if(ina!=inb&&std::abs(da-db)>1e-14){double t=std::clamp(da/(da-db),0.0,1.0);out.push_back({a.q+t*(b.q-a.q),a.bary+t*(b.bary-a.bary)});}}
        return out;};
    auto clipPlane=[&](const std::vector<Clip>& face,const PanelPlane& plane,const Triangle& triangle,bool above){std::vector<Clip> out;if(face.empty())return out;
        auto distance=[&](const Clip& q){V point=V::Zero();for(int j=0;j<3;++j)point+=q.bary[j]*p[triangle.vertices[j]];return (point-plane.point).dot(plane.normal);};
        bool coplanar=true;for(const auto& vertex:face)if(std::abs(distance(vertex))>1e-9){coplanar=false;break;}
        if(coplanar)return above?face:out;
        for(size_t i=0;i<face.size();++i){const auto& a=face[i];const auto& b=face[(i+1)%face.size()];double da=distance(a),db=distance(b);bool ina=above?da>=-1e-10:da<=1e-10,inb=above?db>=-1e-10:db<=1e-10;
            if(ina)out.push_back(a);if(ina!=inb&&std::abs(da-db)>1e-14){double t=std::clamp(da/(da-db),0.0,1.0);out.push_back({a.q+t*(b.q-a.q),a.bary+t*(b.bary-a.bary)});}}
        return out;
    };
    std::map<std::pair<int,int>,UV> nodeUV;
    std::vector<Polygon> output;std::vector<int> bands,columns,outputCharts;std::vector<double> areas;
    for(int ti=0;ti<int(triangles.size());++ti){const auto& triangle=triangles[ti];int chartId=chartIds[triangle.polygon];const auto& chart=charts.at(chartId);
        size_t outputStart=output.size();
        std::vector<Clip> original;UV low=UV::Constant(std::numeric_limits<double>::infinity()),high=-low;
        for(int j=0;j<3;++j){UV q=chart.uv.at(triangle.vertices[j]);const auto& frame=sources[triangle.source].panelFrame;
            if(!frame.empty()){const auto& point=p[triangle.vertices[j]];q=UV(point.dot(V(frame[0],frame[1],frame[2]))-frame[3],point.dot(V(frame[4],frame[5],frame[6]))-frame[7]);}
            if(s.reverse&&!chart.flow)q.x()=-q.x();for(int d=chart.layout?2:(chart.flow?1:0);d<2;++d){double size=d?s.width:s.spacing,offset=d?0:s.offset,level=offset+std::round((q[d]-offset)/size)*size;if(std::abs(q[d]-level)<1e-6*std::min(1.0,size))q[d]=level;}original.push_back({q,V::Unit(j)});low=low.cwiseMin(q);high=high.cwiseMax(q);}
        UV size(s.spacing,s.width),offset(s.offset,0);int first[2],last[2];for(int d=0;d<2;++d){double a=std::floor((low[d]-offset[d])/size[d]),b=std::floor((high[d]-offset[d])/size[d]);if(std::abs(a)>1000000||std::abs(b)>1000000||b-a>1000)throw std::invalid_argument("Panel spacing or width is too small.");first[d]=int(a);last[d]=int(b);}
        if((last[0]-first[0]+1.0)*(last[1]-first[1]+1.0)>10000)throw std::invalid_argument("Too many panel cells per face.");
        auto emit=[&](const std::vector<Clip>& cell,int u,int v){if(cell.size()<3)return;Polygon polygon;polygon.source=triangle.source;
                for(const auto& point:cell){int id=nodeId(point,ti);if(chart.flow)nodeUV[{chartId,id}]=point.q;if(polygon.vertices.empty()||polygon.vertices.back()!=id)polygon.vertices.push_back(id);}if(polygon.vertices.size()>1&&polygon.vertices.front()==polygon.vertices.back())polygon.vertices.pop_back();if(polygon.vertices.size()<3)return;
                double area=0;V origin=nodes[polygon.vertices[0]].p;for(size_t j=1;j+1<polygon.vertices.size();++j)area+=(nodes[polygon.vertices[j]].p-origin).cross(nodes[polygon.vertices[j+1]].p-origin).norm()/2;
                if(area<1e-12)return;output.push_back(std::move(polygon));bands.push_back(u);columns.push_back(v);outputCharts.push_back(chartId);areas.push_back(area);
                if(output.size()>2000000||nodes.size()>2000000)throw std::invalid_argument("Panel output exceeds two million faces/vertices.");
        };
        if(chart.layout){
            for(int v=0;v<int(chart.loopPlanes.size());++v){auto remaining=clip(clip(original,1,chart.loopCuts[v],true),1,chart.loopCuts[v+1],false);if(remaining.size()<3)continue;
                const auto& planes=chart.loopPlanes[v];for(int u=0;u<int(planes.size())&&remaining.size()>=3;++u){emit(clipPlane(remaining,planes[u],triangle,false),u,v);remaining=clipPlane(remaining,planes[u],triangle,true);}emit(remaining,int(planes.size()),v);
            }
        } else if(!chart.flow) {
            for(int u=first[0];u<=last[0];++u){auto strip=clip(clip(original,0,s.offset+u*s.spacing,true),0,s.offset+(u+1)*s.spacing,false);if(strip.size()<3)continue;
                for(int v=first[1];v<=last[1];++v)emit(clip(clip(strip,1,v*s.width,true),1,(v+1)*s.width,false),u,v);
            }
        } else {
            for(int v=first[1];v<=last[1];++v){auto strip=clip(clip(original,1,v*s.width,true),1,(v+1)*s.width,false);if(strip.size()<3)continue;
                std::vector<double> cuts;for(auto cut:chart.flowCuts)if(cut.second<(v+1)*s.width-1e-8||cut.second==0)cuts.push_back(cut.first);
                std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end(),[](double a,double b){return std::abs(a-b)<1e-8;}),cuts.end());
                for(int u=0;u+1<int(cuts.size());++u)if(cuts[u]<=high.x()+1e-9&&cuts[u+1]>=low.x()-1e-9)emit(clip(clip(strip,0,cuts[u],true),0,cuts[u+1],false),u,v);
            }
        }
        if(chart.flow||chart.layout){
            // An adaptive split starts on a Y joint. Insert that junction into both
            // neighbouring cells, including joints lying inside an input triangle.
            std::set<int> triangleNodes;for(size_t f=outputStart;f<output.size();++f)triangleNodes.insert(output[f].vertices.begin(),output[f].vertices.end());
            for(size_t f=outputStart;f<output.size();++f){auto& loop=output[f].vertices;std::vector<int> refined;
                for(size_t j=0;j<loop.size();++j){int a=loop[j],b=loop[(j+1)%loop.size()];refined.push_back(a);V direction=nodes[b].p-nodes[a].p;double length2=direction.squaredNorm();if(length2<1e-20)continue;std::vector<std::pair<double,int>> between;
                    for(int n:triangleNodes)if(n!=a&&n!=b){double t=(nodes[n].p-nodes[a].p).dot(direction)/length2;if(t>1e-8&&t<1-1e-8&&(nodes[n].p-nodes[a].p-t*direction).norm()<1e-8)between.push_back({t,n});}
                    std::sort(between.begin(),between.end());for(auto n:between)refined.push_back(n.second);
                }loop=std::move(refined);
            }
        }
    }
    // Charts meeting at creases may have different subdivisions. Insert the union of
    // seam cuts on both sides, avoiding T-junctions without moving surface geometry.
    auto supportEdge=[&](int a,int b)->Edge{std::set<int> support={nodes[a].a,nodes[a].b,nodes[b].a,nodes[b].b};if(support.count(-1)||support.size()!=2)return {-1,-1};return {*support.begin(),*support.rbegin()};};
    for(auto& polygon:output){std::vector<int> refined;const auto& f=polygon.vertices;for(size_t i=0;i<f.size();++i){int a=f[i],b=f[(i+1)%f.size()];refined.push_back(a);auto e=supportEdge(a,b);auto cuts=edgeCuts.find(e);if(cuts==edgeCuts.end())continue;
            auto parameter=[&](int id){return nodes[id].a==nodes[id].b?(id==e.first?0.0:1.0):nodes[id].t;};double ta=parameter(a),tb=parameter(b);std::vector<std::pair<double,int>> between;
            for(auto item:cuts->second){int id=item.second;double t=nodes[id].t;if(id!=a&&id!=b&&t>std::min(ta,tb)+1e-9&&t<std::max(ta,tb)-1e-9)between.push_back({t,id});}if(tb<ta)std::reverse(between.begin(),between.end());for(auto item:between)refined.push_back(item.second);
        }polygon.vertices=std::move(refined);}
    auto onCrease=[&](int a,int b){std::set<int> support={nodes[a].a,nodes[a].b,nodes[b].a,nodes[b].b};return support.size()==2&&creases.count(edge(*support.begin(),*support.rbegin()));};
    auto onDivider=[&](int a,int b){auto e=supportEdge(a,b);return dividers.count(e)>0;};
    auto flowCoordinate=[&](int chartId,int id)->UV {
        auto found=nodeUV.find({chartId,id});if(found!=nodeUV.end())return found->second;
        const auto& chart=charts[chartId];const auto& n=nodes[id];
        if(n.a>=0&&chart.uv.count(n.a)&&chart.uv.count(n.b))return (1-n.t)*chart.uv.at(n.a)+n.t*chart.uv.at(n.b);
        return UV(-1e30,-1e30);
    };
    auto redRail=[&](int chartId,int a,int b){const auto& chart=charts[chartId];if(!chart.flow)return false;
        auto qa=flowCoordinate(chartId,a),qb=flowCoordinate(chartId,b);
        if(std::abs(qa.x()-qb.x())>1e-7)return false;
        auto it=chart.redCuts.lower_bound(qa.x()-1e-7);return it!=chart.redCuts.end()&&std::abs(*it-qa.x())<1e-7;
    };
    auto beyondBend=[&](int chartId,int id){const auto& chart=charts[chartId];if(chart.bendDepth.empty())return false;auto q=flowCoordinate(chartId,id);
        auto b=chart.bendDepth.lower_bound(int(std::ceil(q.x())));double depth=1e30;
        if(b==chart.bendDepth.end())depth=chart.bendDepth.rbegin()->second;
        else if(b==chart.bendDepth.begin()||std::abs(b->first-q.x())<1e-7)depth=b->second;
        else {auto a=std::prev(b);double t=(q.x()-a->first)/(b->first-a->first);depth=(1-t)*a->second+t*b->second;}
        return q.y()>=depth-1e-7;
    };
    std::map<Edge,std::vector<int>> outputEdges;
    for(int i=0;i<int(output.size());++i){const auto& f=output[i].vertices;for(size_t j=0;j<f.size();++j)outputEdges[edge(f[j],f[(j+1)%f.size()])].push_back(i);}
    DSU panels(int(output.size()));
    for(const auto& item:outputEdges) {
        const auto& f=item.second;if(f.size()>2){std::string detail;for(int face:f)detail+=" source="+std::to_string(output[face].source)+" band="+std::to_string(bands[face])+" column="+std::to_string(columns[face])+" area="+std::to_string(areas[face]);throw std::invalid_argument("Non-manifold clipped output: edge length="+std::to_string((nodes[item.first.first].p-nodes[item.first.second].p).norm())+detail);}
        if(f.size()==2&&bands[f[0]]==bands[f[1]]&&columns[f[0]]==columns[f[1]]&&outputCharts[f[0]]==outputCharts[f[1]]&&output[f[0]].source==output[f[1]].source&&!onCrease(item.first.first,item.first.second)&&!onDivider(item.first.first,item.first.second))panels.join(f[0],f[1]);
    }
    // Stop Y joints one flow cell before a curved divider. Merge the two panel
    // regions sharing that whole joint, rather than leaving a free slit or merely
    // hiding its display line. Divider and flow boundaries remain panel edges.
    if(!dividers.empty()) {
        std::set<int> dividerVertices;
        for(const auto& item:outputEdges)if(onDivider(item.first.first,item.first.second)){dividerVertices.insert(item.first.first);dividerVertices.insert(item.first.second);}
        std::set<int> redVertices;
        for(const auto& item:outputEdges){const auto& f=item.second;auto [a,b]=item.first;
            if(f.size()==1||onCrease(a,b)||onDivider(a,b)||redRail(outputCharts[f[0]],a,b)||
               (f.size()==2&&(outputCharts[f[0]]!=outputCharts[f[1]]||output[f[0]].source!=output[f[1]].source))){redVertices.insert(a);redVertices.insert(b);}
        }
        std::map<Edge,int> groups;std::vector<Edge> pairs;std::vector<std::pair<Edge,int>> yEdges;
        std::set<int> remove;std::map<int,double> lengths;
        for(const auto& item:outputEdges){const auto& f=item.second;auto [a,b]=item.first;
            if(f.size()!=2||columns[f[0]]==columns[f[1]]||outputCharts[f[0]]!=outputCharts[f[1]]||output[f[0]].source!=output[f[1]].source||onCrease(a,b)||onDivider(a,b))continue;
            auto pair=edge(panels.root(f[0]),panels.root(f[1]));auto group=groups.emplace(pair,int(groups.size()));int id=group.first->second;
            if(group.second)pairs.push_back(pair);yEdges.push_back({item.first,id});
            if(!charts[outputCharts[f[0]]].bendDepth.empty())lengths[id]+=(nodes[a].p-nodes[b].p).norm();
            if(dividerVertices.count(a)||dividerVertices.count(b)||beyondBend(outputCharts[f[0]],a)||beyondBend(outputCharts[f[0]],b))remove.insert(id);
        }
        // Continue cleanup through blue flow boundaries, stopping at red rails.
        // Thus the last cyan segment ends on red rather than on an added blue path.
        DSU spans(int(pairs.size()));std::map<int,int> atVertex;
        for(auto item:yEdges)for(int v:{item.first.first,item.first.second})if(!redVertices.count(v)){
            auto found=atVertex.emplace(v,item.second);if(!found.second)spans.join(item.second,found.first->second);
        }
        std::map<int,double> spanLength;for(auto item:lengths)spanLength[spans.root(item.first)]+=item.second;
        std::set<int> removeSpans;for(int id:remove)removeSpans.insert(spans.root(id));
        for(auto item:spanLength)if(item.second<0.1*std::min(s.spacing,s.width))removeSpans.insert(item.first);
        for(int id=0;id<int(pairs.size());++id)if(removeSpans.count(spans.root(id)))panels.join(pairs[id].first,pairs[id].second);
    }
    std::vector<zPanelInfo> info;std::map<int,int> panelMap;zIntArray panelIds,sourceIds,counts,connects,joints,kinds;
    for(int i=0;i<int(output.size());++i) {
        auto found=panelMap.emplace(panels.root(i),int(info.size()));int id=found.first->second;
        if(found.second)info.push_back({output[i].source,bands[i],0,columns[i],outputCharts[i]});info[id].area+=areas[i];
        panelIds.push_back(id);sourceIds.push_back(output[i].source);counts.push_back(int(output[i].vertices.size()));
        connects.insert(connects.end(),output[i].vertices.begin(),output[i].vertices.end());
    }
    for(const auto& item:outputEdges) {
        auto [a,b]=item.first;const auto& f=item.second;
        if(f.size()==2&&panelIds[f[0]]==panelIds[f[1]])continue;
        joints.push_back(a);joints.push_back(b);
        kinds.push_back(f.size()==1?0:onDivider(a,b)?4:outputCharts[f[0]]!=outputCharts[f[1]]||onCrease(a,b)||redRail(outputCharts[f[0]],a,b)||output[f[0]].source!=output[f[1]].source?2:columns[f[0]]!=columns[f[1]]?3:1);
    }
    zPointArray vertices;for(const auto& n:nodes)vertices.emplace_back(float(n.p.x()),float(n.p.y()),float(n.p.z()));
    zObjectMesh mesh;zFnMesh(mesh).create(vertices,counts,connects);
    m_guides=std::move(guides);m_mesh=std::move(mesh);m_facePanels=std::move(panelIds);m_faceSources=std::move(sourceIds);m_jointEdges=std::move(joints);m_jointKinds=std::move(kinds);m_panels=std::move(info);
    m_components=int(groups.size());m_weldedEdges=welded;m_axis=0;m_smallPanels=0;for(const auto& panel:m_panels)if(panel.area<s.smallPanelArea)++m_smallPanels;
}
}
