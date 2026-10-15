// C-Tubes construction adapted from EPFL-LGG/C-tubes (MIT).
// See docs/third-party/C-tubes-LICENSE.txt for attribution.
#include <zspace/zToolsets/fabrication/zTsCTubes.h>
#include <Eigen/Dense>
#include "../../../third_party/LBFGSpp/include/LBFGSB.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <set>
#include <stdexcept>

namespace zSpace
{
namespace
{
using V = Eigen::Vector3d;
V vec(const zPoint& p) { return V(p.x, p.y, p.z); }
zPoint point(const V& p) { return zPoint(float(p.x()), float(p.y()), float(p.z())); }
V unit(const V& v)
{
    if (!v.allFinite() || v.norm() < 1e-10) throw std::invalid_argument("Degenerate edge or reversing tangent.");
    return v.normalized();
}
V project(const V& v, const V& n, const V& direction)
{
    const double d = n.dot(direction);
    if (std::abs(d) < 1e-10 * direction.norm()) throw std::invalid_argument("Tube projection is parallel to its section plane.");
    if (direction.norm() < 1e-12) throw std::invalid_argument("Collapsed tube projection direction.");
    return v - direction * (n.dot(v) / d);
}
// An open chain has a plane when bent. A straight line uses the path direction
// projected perpendicular to the line; its roll is otherwise underdetermined.
V profileNormal(const std::vector<V>& points, const V& preferred) {
    // Use all vertices, independent of graph traversal or branch order.
    V normal=V::Zero(),direction=V::Zero();double best=0;
    for(size_t i=1;i<points.size();++i)if((points[i]-points[0]).squaredNorm()>best){direction=points[i]-points[0];best=direction.squaredNorm();}
    best=0;
    for(size_t j=1;j<points.size();++j) {
        V cross=direction.cross(points[j]-points[0]);
        if(cross.squaredNorm()>best){best=cross.squaredNorm();normal=cross;}
    }
    if(normal.norm()>1e-10) return unit(normal);
    V x=unit(points[1]-points[0]);normal=preferred-x*x.dot(preferred);
    if(normal.norm()<1e-10) {V axis=std::abs(x.z())<.9?V::UnitZ():V::UnitY();normal=axis-x*x.dot(axis);}
    return unit(normal);
}
zPointArray orderedGraph(zObjectGraph& graph, bool& closed)
{
    zFnGraph fn(graph); zPointArray p; zIntArray e;
    fn.getVertexPositions(p); fn.getEdgeData(e);
    if (p.size() < 2 || e.empty() || e.size() % 2) throw std::invalid_argument("Assign a connected polyline graph.");
    std::vector<std::vector<int>> adjacency(p.size());
    std::set<std::pair<int,int>> edges;
    for (size_t i = 0; i < e.size(); i += 2)
    {
        int a=e[i], b=e[i+1];
        if (a<0 || b<0 || a>=int(p.size()) || b>=int(p.size()) || a==b || !edges.insert(std::minmax(a,b)).second)
            throw std::invalid_argument("Invalid or duplicate graph edge.");
        if ((vec(p[a])-vec(p[b])).norm()<1e-10) throw std::invalid_argument("Zero length graph edge.");
        adjacency[a].push_back(b); adjacency[b].push_back(a);
    }
    int endpoints=0, start=0;
    for (int i=0; i<int(p.size()); ++i)
    {
        if (!vec(p[i]).allFinite() || adjacency[i].empty() || adjacency[i].size()>2)
            throw std::invalid_argument("Graph must be finite, connected and unbranched.");
        if (adjacency[i].size()==1) { if (!endpoints) start=i; ++endpoints; }
    }
    if (endpoints!=0 && endpoints!=2) throw std::invalid_argument("Graph is not a single chain or loop.");
    closed=endpoints==0;
    zPointArray out; std::vector<bool> visited(p.size(),false);
    int previous=-1, current=start;
    while (!visited[current])
    {
        visited[current]=true; out.push_back(p[current]);
        int next=-1;
        for (int candidate:adjacency[current]) if(candidate!=previous) { next=candidate; break; }
        if (next<0) break;
        previous=current; current=next;
    }
    if (out.size()!=p.size()) throw std::invalid_argument("Disconnected graph.");
    if (closed) out.push_back(out.front());
    return out;
}
void mesh(zObjectMesh& out, zPointArray p, zIntArray counts, zIntArray connects);
// Preserve branch vertex IDs. Chains retain the historical ordered layout.
zPointArray sectionGraph(zObjectGraph& graph,zIntArray& edges,bool& closed,bool& branched) {
    zFnGraph fn(graph);zPointArray points;fn.getVertexPositions(points);fn.getEdgeData(edges);
    if(points.size()<2 || edges.empty() || edges.size()%2) throw std::invalid_argument("Supply a connected section graph.");
    std::vector<std::vector<int>> adjacency(points.size());std::set<std::pair<int,int>> unique;
    for(size_t k=0;k<edges.size();k+=2){int a=edges[k],b=edges[k+1];
        if(a<0||b<0||a>=int(points.size())||b>=int(points.size())||a==b||!unique.insert(std::minmax(a,b)).second)
            throw std::invalid_argument("Invalid or duplicate section edge.");
        if((vec(points[a])-vec(points[b])).norm()<1e-10)throw std::invalid_argument("Zero length section edge.");
        adjacency[a].push_back(b);adjacency[b].push_back(a);
    }
    branched=false;for(size_t i=0;i<points.size();++i){if(!vec(points[i]).allFinite()||adjacency[i].empty())throw std::invalid_argument("Invalid section vertex.");branched=branched||adjacency[i].size()>2;}
    std::set<int> seen;std::vector<int> queue={0};while(!queue.empty()){int i=queue.back();queue.pop_back();if(!seen.insert(i).second)continue;for(int j:adjacency[i])queue.push_back(j);}
    if(seen.size()!=points.size())throw std::invalid_argument("Disconnected section graph.");
    if(branched){closed=false;return points;}
    points=orderedGraph(graph,closed);if(closed)points.pop_back();edges.clear();
    for(int i=0;i<int(points.size())-1;++i)edges.insert(edges.end(),{i,i+1});if(closed)edges.insert(edges.end(),{int(points.size())-1,0});return points;
}
std::vector<zIntArray> sectionMappings(int n,const zIntArray& a,const zIntArray& b) {
    std::vector<std::vector<bool>> aa(n,std::vector<bool>(n)),bb=aa;std::vector<int> da(n),db(n);
    for(size_t k=0;k<a.size();k+=2){aa[a[k]][a[k+1]]=aa[a[k+1]][a[k]]=true;++da[a[k]];++da[a[k+1]];}
    for(size_t k=0;k<b.size();k+=2){bb[b[k]][b[k+1]]=bb[b[k+1]][b[k]]=true;++db[b[k]];++db[b[k+1]];}
    zIntArray order(n),mapping(n,-1);for(int i=0;i<n;++i)order[i]=i;
    std::stable_sort(order.begin(),order.end(),[&](int i,int j){return da[i]>da[j];});
    std::vector<bool> used(n);std::vector<zIntArray> results;int trials=0;
    std::function<void(int)> visit=[&](int depth){
        if(++trials>100000 || results.size()>4096)throw std::invalid_argument("Section topology has too many ambiguous correspondences; simplify branches.");
        if(depth==n){results.push_back(mapping);return;}int i=order[depth];
        for(int j=0;j<n;++j)if(!used[j]&&da[i]==db[j]){bool valid=true;for(int k=0;k<depth;++k){int v=order[k];if(aa[i][v]!=bb[j][mapping[v]])valid=false;}
            if(valid){mapping[i]=j;used[j]=true;visit(depth+1);used[j]=false;mapping[i]=-1;}}
    };visit(0);if(results.empty())throw std::invalid_argument("Target section graphs must have matching connectivity.");return results;
}
// A branched sweep is nonmanifold at its longitudinal junction. zObjectMesh
// uses manifold halfedges, so render coincident strip seams with separate IDs.
void sectionMesh(zObjectMesh& out,const zPointArray& points,const zIntArray& edges,int stations,bool branched) {
    int n=int(points.size())/stations;zPointArray vertices=branched?zPointArray{}:points;zIntArray counts,connects;
    for(size_t k=0;k<edges.size();k+=2){int a=edges[k],b=edges[k+1],base=int(vertices.size());
        if(branched)for(int j=0;j<stations;++j){vertices.push_back(points[j*n+a]);vertices.push_back(points[j*n+b]);}
        for(int j=0;j<stations-1;++j){counts.push_back(4);if(branched)connects.insert(connects.end(),{base+2*j,base+2*j+2,base+2*j+3,base+2*j+1});else connects.insert(connects.end(),{j*n+a,(j+1)*n+a,(j+1)*n+b,j*n+b});}
    }mesh(out,vertices,counts,connects);
}
std::vector<V> sample(const zPointArray& points, bool closed, const zCTubeSettings& settings)
{
    std::vector<V> p; for (const auto& v:points) p.push_back(vec(v));
    if (!settings.spline) return p;
    const int n=int(p.size()), unique=closed?n-1:n;
    if(unique>512) throw std::invalid_argument("Spline mode supports up to 512 control points; use polyline mode for denser graphs.");
    if (closed && unique<3) throw std::invalid_argument("A closed directrix needs at least three vertices.");
    Eigen::MatrixXd A=Eigen::MatrixXd::Zero(unique,unique), rhs=Eigen::MatrixXd::Zero(unique,3);
    for (int i=0; i<unique; ++i)
    {
        if (!closed && (i==0 || i==unique-1)) { A(i,i)=1; continue; }
        int prev=(i+unique-1)%unique, next=(i+1)%unique;
        A(i,i)=4; A(i,prev)+=1; A(i,next)+=1;
        rhs.row(i)=(6.0*(p[next]-2.0*p[i]+p[prev])).transpose();
    }
    // Natural boundary rows make the open system non-symmetric.
    Eigen::MatrixXd second=A.partialPivLu().solve(rhs);
    std::vector<V> result; result.reserve(settings.samples);
    for (int i=0; i<settings.samples; ++i)
    {
        double t=double(i)*(n-1)/(settings.samples-1);
        int k=std::min(int(t),n-2); double b=t-k, a=1-b;
        V s0=second.row(k%unique).transpose(), s1=second.row((k+1)%unique).transpose();
        result.push_back(a*p[k]+b*p[k+1]+((a*a*a-a)*s0+(b*b*b-b)*s1)/6.0);
    }
    return result;
}
void mesh(zObjectMesh& out, zPointArray p, zIntArray counts, zIntArray connects)
{
    zFnMesh fn(out); fn.create(p,counts,connects); fn.computeMeshNormals();
}
}
void zTsCTubes::clearResults()
{
    m_fit=zCTubeFitReport();
    m_frames.clear(); m_stripAreas.clear();
    m_vertices.clear(); m_stations=0; m_strips.clear(); m_unrolled.clear();
    m_tube=zObjectMesh(); m_sampledDirectrix=zObjectGraph();
}
void zTsCTubes::setDirectrix(zObjectGraph& graph)
{
    bool closed; auto p=orderedGraph(graph,closed); m_curve=std::move(p); m_closed=closed; clearResults();
}
void zTsCTubes::setGeneratrix(zObjectGraph& graph)
{
    bool closed,branched;zIntArray edges;auto p=sectionGraph(graph,edges,closed,branched);
    if(p.size()<2) throw std::invalid_argument("Generatrix needs at least two distinct vertices.");
    V center=V::Zero(); for(const auto& v:p) center+=vec(v); center/=double(p.size());
    std::vector<V> world;for(auto v:p)world.push_back(vec(v));
    V preferred=m_curve.size()>1?unit(vec(m_curve[1])-vec(m_curve[0])):V::UnitZ();
    V normal=profileNormal(world,preferred); V x=unit(vec(p[1])-vec(p[0])), y=unit(normal.cross(x));
    double extent=0; for(const auto& v:p) extent=std::max(extent,(vec(v)-center).norm());
    // zPoint stores floats: account for coordinate quantisation in georeferenced inputs.
    double coordinateMagnitude=0; for(const auto& v:p) coordinateMagnitude=std::max(coordinateMagnitude,vec(v).cwiseAbs().maxCoeff());
    const double planeTolerance=1e-5*std::max(extent,1e-6)+4*std::numeric_limits<float>::epsilon()*coordinateMagnitude;
    zPointArray profile;
    for(const auto& v:p)
    {
        V d=vec(v)-center;
        if(std::abs(d.dot(normal))>planeTolerance) throw std::invalid_argument("Generatrix must be planar.");
        profile.emplace_back(float(d.dot(x)),float(d.dot(y)),0.0f);
    }
    m_worldProfile=p; m_targets.clear();
    m_origins.resize(1); m_hasOrigin.assign(1,false);
    m_profileX=point(x); m_profileY=point(y); m_profileNormal=point(normal);
    m_profile=std::move(profile); m_profileClosed=closed;m_profileBranched=branched;m_profileEdges=std::move(edges);m_targetMappings.clear(); clearResults();
}
void zTsCTubes::setGeneratrices(zObjectGraphArray& graphs)
{
    if(graphs.empty() || graphs.size()>16) throw std::invalid_argument("Supply 1 to 16 generatrices.");
    std::vector<zPointArray> targets;bool targetClosed=true;std::vector<zIntArray> targetEdges;bool branched=false;
    for(auto& graph:graphs) {
        zTsCTubes validator; validator.setGeneratrix(graph);
        if(!targets.empty() && validator.m_profileClosed!=targetClosed) throw std::invalid_argument("Target profiles must have matching open or closed topology.");
        if(targets.empty()){targetClosed=validator.m_profileClosed;branched=validator.m_profileBranched;}
        if(branched!=validator.m_profileBranched)throw std::invalid_argument("Target sections must have matching connectivity.");
        targetEdges.push_back(validator.m_profileEdges);
        if(!targets.empty() && validator.m_worldProfile.size()!=targets[0].size())
            throw std::invalid_argument("Target generatrices need matching vertex counts; use matching curve sampling or polygon topology.");
        targets.push_back(validator.m_worldProfile);
    }
    std::vector<std::vector<zIntArray>> mappings;
    if(branched)for(auto& e:targetEdges)mappings.push_back(sectionMappings(int(targets[0].size()),targetEdges[0],e));
    setGeneratrix(graphs.front()); m_targets=std::move(targets);m_targetMappings=std::move(mappings);
    m_origins.resize(graphs.size()); m_hasOrigin.assign(graphs.size(),false);
}
void zTsCTubes::setGeneratrixOrigin(int index, const zPoint& vertex)
{
    if(index<0 || index>=int(m_hasOrigin.size())) throw std::invalid_argument("Invalid generatrix origin index.");
    const auto& points=m_targets.empty()?m_worldProfile:m_targets[index];
    bool found=false; for(auto p:points) if((vec(p)-vec(vertex)).norm()<1e-7) found=true;
    if(!found) throw std::invalid_argument("Generatrix origin must be one of its vertices.");
    m_origins[index]=vertex; m_hasOrigin[index]=true; clearResults();
}
void zTsCTubes::setSettings(const zCTubeSettings& s)
{
    auto validWeight=[](double w){return std::isfinite(w) && w>=0 && w<=1e8;};
    for(double w:{s.targetSectionWeight,s.targetPlaneWeight,s.apexSmoothnessWeight,s.planeSmoothnessWeight})
        if(!validWeight(w)) throw std::invalid_argument("Objective weights must be finite and in [0, 1e8].");
    for(const auto& weights:{s.sectionWeights,s.planeWeights}) {
        if(!weights.empty() && weights.size()!=m_hasOrigin.size()) throw std::invalid_argument("Per-target weights must match the generatrix count.");
        for(double w:weights) if(!validWeight(w)) throw std::invalid_argument("Invalid per-target weight.");
    }
    if(s.fitIterations<1 || s.fitIterations>2000 || !std::isfinite(s.fitTolerance) || s.fitTolerance<=0) throw std::invalid_argument("Invalid fit settings.");
    if(s.samples<2 || s.samples>4096 || !std::isfinite(s.scale) || s.scale<=0 ||
       !std::isfinite(s.twistRadians) || !std::isfinite(s.apexScale) || s.apexScale<=0)
        throw std::invalid_argument("Invalid tube settings (samples 2..4096; scales positive).");
    m_settings=s; clearResults();
}
void zTsCTubes::computeTube(zCTubeMode mode)
{
    clearResults();
    if(m_curve.empty() || m_profile.empty()) throw std::invalid_argument("Assign both directrix and generatrix graphs.");
    if(mode!=zCTubeMode::Generic && mode!=zCTubeMode::ParallelTransport) throw std::invalid_argument("Unknown tube mode.");
    if(m_targets.size()>1) {
        if(mode!=zCTubeMode::Generic) throw std::invalid_argument("Target-section optimisation uses Generic C-Tube; parallel transport has no apex-locating function.");
        computeFittedTube(); return;
    }
    auto X=sample(m_curve,m_closed,m_settings); int M=int(X.size()), N=int(m_profile.size());
    if(size_t(M)*N>2000000) throw std::invalid_argument("Tube exceeds two million vertices.");
    std::vector<V> t,n(M),v(size_t(M)*N);
    for(int i=1;i<M;++i) t.push_back(unit(X[i]-X[i-1]));
    n.front()=m_closed?unit(t.front()+t.back()):t.front(); n.back()=m_closed?n.front():t.back();
    for(int i=1;i<M-1;++i) n[i]=unit(t[i-1]+t[i]);
    V axis=std::abs(t[0].z())<0.9?V::UnitZ():V::UnitY();
    V e1=unit(t[0].cross(axis)), e2=t[0].cross(e1);
    if(m_settings.preserveOrientation)
    {
        e1=vec(m_profileX); e2=vec(m_profileY); n[0]=vec(m_profileNormal);
        if(std::abs(n[0].dot(t[0]))<1e-8) throw std::invalid_argument("Generatrix plane is parallel to the initial path; rotate the profile or disable Preserve orientation.");
    }
    double c=std::cos(m_settings.twistRadians), s=std::sin(m_settings.twistRadians);
    std::vector<V> offsets;
    V localOrigin=V::Zero();
    if(m_hasOrigin[0]) {
        V center=V::Zero(); for(auto p:m_worldProfile) center+=vec(p); center/=N;
        V d=vec(m_origins[0])-center; localOrigin=V(d.dot(vec(m_profileX)),d.dot(vec(m_profileY)),0);
    }
    for(int k=0;k<N;++k)
    {
        V p=(vec(m_profile[k])-localOrigin)*m_settings.scale;
        V d=e1*(c*p.x()-s*p.y())+e2*(s*p.x()+c*p.y()); offsets.push_back(d);
        v[k]=X[0]+project(d,n[0],t[0]);
    }
    for(int j=1;j<M;++j)
    {
        if(mode==zCTubeMode::ParallelTransport)
        {
            if(j>1)
            {
                V a=t[j-2], b=t[j-1], cross=a.cross(b); double cosine=a.dot(b);
                if(cosine<-1+1e-10) throw std::invalid_argument("Parallel transport cannot cross a reversing tangent.");
                for(auto& d:offsets) d=d+cross.cross(d)+cross.cross(cross.cross(d))/(1+cosine);
            }
            for(int k=0;k<N;++k) v[j*N+k]=X[j]+project(offsets[k],n[j],t[j-1]);
        }
        else
        {
            V op=X[j]+project(X[j-1]-X[j],n[j],t[j-1]);
            for(int k=0;k<N;++k)
            {
                V prev=v[(j-1)*N+k], prime=m_settings.apexScale*(prev-X[j-1])+op;
                v[j*N+k]=op+project(prime-op,n[j],prime-prev);
            }
        }
    }
    std::vector<zCTubeFrame> frames;
    for(int j=0;j<M;++j)
    {
        V x=unit(v[j*N+1]-v[j*N]);
        V normal=unit(n[j]); V y=unit(normal.cross(x));
        frames.push_back({point(X[j]),point(x),point(y),point(normal)});
    }
    zPointArray points; for(const auto& p:v)
    { if(!p.allFinite() || p.cwiseAbs().maxCoeff()>1e20) throw std::runtime_error("Tube construction diverged."); points.push_back(point(p)); }
    sectionMesh(m_tube,points,m_profileEdges,M,m_profileBranched);
    zPointArray curve; zIntArray edges;
    for(int j=0;j<M;++j) { curve.push_back(point(X[j])); if(j) edges.insert(edges.end(),{j-1,j}); }
    zFnGraph(m_sampledDirectrix).create(curve,edges);
    m_vertices=std::move(points); m_stations=M; m_frames=std::move(frames);
}
void zTsCTubes::computeTubeWithFrames(const std::vector<zCTubeFrame>& supplied,
    const zPointArray& flangeTargets, int flangeVertex)
{
    clearResults();
    if(m_closed || m_worldProfile.empty() || m_targets.size()>1 || supplied.size()!=m_curve.size() || supplied.size()<2)
        throw std::invalid_argument("Fixed frames require one seed profile and one frame per open directrix vertex.");
    const int M=int(supplied.size()),N=int(m_worldProfile.size());
    if(size_t(M)*N>2000000)throw std::invalid_argument("Tube exceeds two million vertices.");
    std::vector<V> X,n,v(size_t(M)*N);
    for(int j=0;j<M;++j) {
        V origin=vec(supplied[j].origin),x=unit(vec(supplied[j].x)),y=unit(vec(supplied[j].y));
        if((origin-vec(m_curve[j])).norm()>1e-5*std::max(1.0,origin.norm()))
            throw std::invalid_argument("Fixed frame origin must be its directrix vertex.");
        if(std::abs(x.dot(y))>1e-5)throw std::invalid_argument("Fixed frame axes must be perpendicular.");
        V normal=unit(x.cross(y));
        if(normal.dot(unit(vec(supplied[j].normal)))<1-1e-5)throw std::invalid_argument("Invalid fixed frame normal.");
        X.push_back(origin);n.push_back(normal);
    }
    for(int k=0;k<N;++k) {
        v[k]=vec(m_worldProfile[k]);
        if(std::abs(n[0].dot(v[k]-X[0]))>1e-5*std::max(1.0,(v[k]-X[0]).norm()))
            throw std::invalid_argument("Seed profile must lie in the first fixed frame plane.");
    }
    for(int j=1;j<M;++j) {
        V direction=unit(X[j]-X[j-1]);
        V op=X[j]+project(X[j-1]-X[j],n[j],direction);
        for(int k=0;k<N;++k) {
            V prev=v[(j-1)*N+k],prime=m_settings.apexScale*(prev-X[j-1])+op;
            v[j*N+k]=op+project(prime-op,n[j],prime-prev);
        }
    }
    if(!flangeTargets.empty()) {
        if(flangeTargets.size()!=supplied.size() || flangeVertex<0 || flangeVertex>=N)
            throw std::invalid_argument("Supply one flange target per frame and a valid flange vertex.");
        int anchor=-1,outer=-1;
        for(int k=0;k<N;++k)if((vec(m_worldProfile[k])-X[0]).norm()<1e-5)anchor=k;
        if(anchor<0 || anchor==flangeVertex)throw std::invalid_argument("Flange closure requires its junction vertex on the directrix.");
        V seed=unit(vec(m_worldProfile[flangeVertex])-X[0]);bool connected=false;
        for(size_t e=0;e<m_profileEdges.size();e+=2) {
            int a=m_profileEdges[e],b=m_profileEdges[e+1],other=a==anchor?b:b==anchor?a:-1;
            if(other==flangeVertex)connected=true;
            if(other>=0 && other!=flangeVertex && seed.dot(unit(vec(m_worldProfile[other])-X[0]))<-1+1e-5)outer=other;
        }
        if(!connected || outer<0)throw std::invalid_argument("Flange closure requires two opposite flange edges at the junction.");
        for(int j=0;j<M;++j) {
            V target=vec(flangeTargets[j]),delta=target-X[j];
            if(!target.allFinite() || delta.norm()<1e-8 || std::abs(n[j].dot(delta))>1e-5*std::max(1.0,delta.norm()))
                throw std::invalid_argument("Flange target must be distinct from the anchor and lie in its frame plane.");
            V x=unit(delta);
            if(x.dot(unit(vec(supplied[j].x)))<1-1e-5)throw std::invalid_argument("Flange target must lie on the inward frame X axis.");
            const double outerLength=(v[j*N+outer]-X[j]).norm();
            v[j*N+anchor]=X[j];v[j*N+flangeVertex]=target;
            // Retain the outward projected length, and a single straight crossbar.
            v[j*N+outer]=X[j]-x*outerLength;
        }
    }
    zPointArray points;for(const auto& p:v) {
        if(!p.allFinite() || p.cwiseAbs().maxCoeff()>1e20)throw std::runtime_error("Fixed-frame construction diverged.");
        points.push_back(point(p));
    }
    sectionMesh(m_tube,points,m_profileEdges,M,m_profileBranched);
    zIntArray edges;for(int j=1;j<M;++j)edges.insert(edges.end(),{j-1,j});
    zFnGraph(m_sampledDirectrix).create(m_curve,edges);
    m_vertices=std::move(points);m_stations=M;m_frames=supplied;
}
void zTsCTubes::computeTubeWithSections(const std::vector<zCTubeFrame>& frames,const zPointArray& sections)
{
    clearResults();const int M=int(frames.size()),N=int(m_worldProfile.size());
    if(M<2 || N<2 || M!=int(m_curve.size()) || sections.size()!=size_t(M)*N || sections.size()>2000000)
        throw std::invalid_argument("Invalid prescribed tube sections.");
    for(int j=0;j<M;++j) {
        V origin=vec(frames[j].origin),normal=unit(vec(frames[j].normal));
        if((origin-vec(m_curve[j])).norm()>1e-5)throw std::invalid_argument("Section origin must match its directrix.");
        for(int k=0;k<N;++k){V p=vec(sections[j*N+k]);if(!p.allFinite() || std::abs(normal.dot(p-origin))>1e-5*std::max(1.0,(p-origin).norm()))throw std::invalid_argument("Prescribed section must lie in its frame plane.");}
    }
    sectionMesh(m_tube,sections,m_profileEdges,M,m_profileBranched);
    zIntArray edges;for(int j=1;j<M;++j)edges.insert(edges.end(),{j-1,j});zFnGraph(m_sampledDirectrix).create(m_curve,edges);
    m_vertices=sections;m_stations=M;m_frames=frames;
}
void zTsCTubes::computeFittedTube()
{
    if(m_closed) throw std::invalid_argument("Target-section fitting currently requires an open directrix.");
    auto path=sample(m_curve,false,m_settings);
    const int N=int(m_targets[0].size()), G=int(m_targets.size());
    if(path.size()>128 || N>128) throw std::invalid_argument("Target fitting supports up to 128 path stations and 128 profile vertices. Reduce sampling.");
    struct Guide { double parameter; V normal; std::vector<V> points; int source; };
    std::vector<Guide> guides;
    for(int g=0;g<G;++g) {
        V center=V::Zero(), normal=V::Zero(); std::vector<V> points;
        for(auto p:m_targets[g]) { points.push_back(vec(p)); center+=vec(p); } center/=N;
        double nearest=std::numeric_limits<double>::infinity();V preferred=V::UnitZ();
        for(size_t j=0;j+1<path.size();++j) {V d=path[j+1]-path[j];double u=std::clamp(d.dot(center-path[j])/d.squaredNorm(),0.0,1.0);double distance=(path[j]+u*d-center).squaredNorm();if(distance<nearest){nearest=distance;preferred=unit(d);}}
        normal=profileNormal(points,preferred);
        double extent=0, magnitude=0;
        for(auto p:points) {extent=std::max(extent,(p-center).norm());magnitude=std::max(magnitude,p.cwiseAbs().maxCoeff());}
        const double planeTolerance=1e-5*std::max(extent,1e-6)+4*std::numeric_limits<float>::epsilon()*magnitude;
        double best=1e100, parameter=-1;
        for(int j=0;j<int(path.size())-1;++j) {
            V d=path[j+1]-path[j]; double denominator=normal.dot(d);
            if(std::abs(denominator)<1e-10*d.norm()) continue;
            double u=normal.dot(center-path[j])/denominator;
            // A path endpoint on a Rhino profile may round just outside its plane.
            double parameterTolerance=std::max(1e-6,planeTolerance/std::abs(denominator));
            if(u<-parameterTolerance || u>1+parameterTolerance) continue;
            u=std::clamp(u,0.0,1.0); double distance=(path[j]+u*d-center).squaredNorm();
            if(distance<best) { best=distance; parameter=j+u; }
        }
        if(parameter<0) throw std::invalid_argument("Each target section plane must intersect the directrix. Extend the path or reposition the section.");
        if(m_hasOrigin[g]) {
            int j=std::min(int(parameter),int(path.size())-2); double u=parameter-j;
            V shift=path[j]*(1-u)+path[j+1]*u-vec(m_origins[g]);
            for(auto& p:points) p+=shift;
        }
        guides.push_back({parameter,normal,points,g});
    }
    std::sort(guides.begin(),guides.end(),[](const Guide& a,const Guide& b){return a.parameter<b.parameter;});
    for(int g=1;g<G;++g) if(guides[g].parameter-guides[g-1].parameter<1e-5)
        throw std::invalid_argument("Two target sections map to the same path station.");
    // Keep exact target intersections, replacing nearly coincident display
    // stations. Tiny endpoint slivers can become shorter than the unroller's
    // global metric tolerance after projection and float storage.
    double pathExtent=0,coordinateMagnitude=0;
    for(auto p:path){pathExtent=std::max(pathExtent,(p-path.front()).norm());coordinateMagnitude=std::max(coordinateMagnitude,p.cwiseAbs().maxCoeff());}
    const double stationTolerance=1e-5*std::max(pathExtent,1e-6)+4*std::numeric_limits<float>::epsilon()*coordinateMagnitude;
    std::vector<double> parameters;
    for(int j=0;j<int(path.size());++j) {
        bool near=false;for(auto& guide:guides){
            int edge=std::min(int(guide.parameter),int(path.size())-2);double u=guide.parameter-edge;
            V intersection=(1-u)*path[edge]+u*path[edge+1];
            if(std::abs(guide.parameter-j)<1 && (intersection-path[j]).norm()<=stationTolerance)near=true;
        }
        if(!near) parameters.push_back(j);
    }
    for(auto& guide:guides) parameters.push_back(guide.parameter);
    std::sort(parameters.begin(),parameters.end());
    const int M=int(parameters.size()); std::vector<V> X,t,n(M);
    for(double parameter:parameters) {
        int j=std::min(int(parameter),int(path.size())-2); double u=parameter-j;
        X.push_back(path[j]*(1-u)+path[j+1]*u);
    }
    for(int j=1;j<M;++j) t.push_back(unit(X[j]-X[j-1]));
    n[0]=t[0]; n[M-1]=t.back(); for(int j=1;j<M-1;++j) n[j]=unit(t[j-1]+t[j]);
    std::vector<int> stations;
    for(auto& guide:guides) {
        int station=int(std::min_element(parameters.begin(),parameters.end(),[&](double a,double b){return std::abs(a-guide.parameter)<std::abs(b-guide.parameter);})-parameters.begin());
        stations.push_back(station);
        if(guide.normal.dot(n[station])<0) guide.normal=-guide.normal;
        n[station]=guide.normal;
    }
    if(m_profileBranched){
        // Keep input-zero topology for the sweep even when another section is first along the path.
        for(auto& guide:guides){const auto original=guide.points;const auto& map=m_targetMappings[guide.source].front();for(int k=0;k<N;++k)guide.points[k]=original[map[k]];}
    }
    int seed=stations[0]; std::vector<int> freeNormals;
    std::vector<V> bx(M),by(M);
    for(int j=0;j<M;++j) {
        if(j!=seed && m_settings.fitPlaneNormals) freeNormals.push_back(j);
        bx[j]=unit(n[j].cross(std::abs(n[j].z())<.9?V::UnitZ():V::UnitY())); by[j]=n[j].cross(bx[j]);
    }
    const int P=M-1+2*int(freeNormals.size());
    Eigen::VectorXd controls=Eigen::VectorXd::Zero(P);
    std::vector<std::vector<V>> targets; for(auto& g:guides) targets.push_back(g.points);
    auto sweep=[&](const Eigen::VectorXd& q,std::vector<V>& vertices,std::vector<V>* outputNormals=nullptr, std::vector<Eigen::MatrixXd>* derivatives=nullptr, std::vector<Eigen::MatrixXd>* normalDerivatives=nullptr) {
        std::vector<V> normals=n;
        if(derivatives) { derivatives->assign(size_t(M)*N,Eigen::MatrixXd::Zero(3,P)); normalDerivatives->assign(M,Eigen::MatrixXd::Zero(3,P)); }
        for(int k=0;k<int(freeNormals.size());++k) {
            int j=freeNormals[k]; V raw=n[j]+bx[j]*q[M-1+2*k]+by[j]*q[M+2*k]; normals[j]=unit(raw);
            if(derivatives) {
                Eigen::Matrix3d dn=(Eigen::Matrix3d::Identity()-normals[j]*normals[j].transpose())/raw.norm();
                (*normalDerivatives)[j].col(M-1+2*k)=dn*bx[j];
                (*normalDerivatives)[j].col(M+2*k)=dn*by[j];
            }
        }
        vertices.resize(size_t(M)*N);
        for(int k=0;k<N;++k) vertices[seed*N+k]=targets[0][k];
        auto step=[&](int from,int to) {
            int edge=std::min(from,to); double a=m_settings.apexScale*std::exp(q[edge]); if(to<from) a=1/a;
            for(int k=0;k<N;++k) {
                V previous=vertices[from*N+k], prime=a*(previous-X[from])+X[to];
                vertices[to*N+k]=X[to]+project(prime-X[to],normals[to],prime-previous);
                if(derivatives) {
                    const auto& dp=(*derivatives)[from*N+k]; const auto& dn=(*normalDerivatives)[to];
                    Eigen::MatrixXd dprime=a*dp;
                    dprime.col(edge)+=(to<from?-a:a)*(previous-X[from]);
                    V v=prime-X[to], direction=prime-previous;
                    double numerator=normals[to].dot(v), denominator=normals[to].dot(direction), ratio=numerator/denominator;
                    Eigen::MatrixXd dd=dprime-dp;
                    Eigen::RowVectorXd dr=((v.transpose()*dn+normals[to].transpose()*dprime)-ratio*(direction.transpose()*dn+normals[to].transpose()*dd))/denominator;
                    (*derivatives)[to*N+k]=dprime-ratio*dd-direction*dr;
                }
                if(!vertices[to*N+k].allFinite() || vertices[to*N+k].norm()>1e12) throw std::runtime_error("Divergent fit trial.");
            }
        };
        for(int j=seed+1;j<M;++j) step(j-1,j);
        for(int j=seed-1;j>=0;--j) step(j+1,j);
        if(outputNormals) *outputNormals=std::move(normals);
    };
    auto evaluate=[&](const Eigen::VectorXd& q,Eigen::VectorXd& residual,double& rms,double& maximum,zCTubeFitReport* detail=nullptr, Eigen::VectorXd* gradient=nullptr) {
        try {
            std::vector<V> vertices,normals; std::vector<Eigen::MatrixXd> dv,dn;
            sweep(q,vertices,&normals,gradient?&dv:nullptr,gradient?&dn:nullptr);
            // Re-evaluate allowed cyclic pairings at every objective call, as upstream.
            for(int g=1;g<G;++g) {
                if(m_profileBranched){
                    double best=std::numeric_limits<double>::infinity();const zIntArray* selected=nullptr;
                    // Seed vertex IDs define the shared topology, independent of input order.
                    const zIntArray& seedMap=m_targetMappings[guides[g].source].front();
                    zIntArray inverse(N);for(int k=0;k<N;++k)inverse[seedMap[k]]=k;
                    for(const auto& map:m_targetMappings[guides[g].source]){double cost=0;for(int k=0;k<N;++k)cost+=(vertices[stations[g]*N+k]-guides[g].points[inverse[map[k]]]).squaredNorm();if(cost<best){best=cost;selected=&map;}}
                    for(int k=0;k<N;++k)targets[g][k]=guides[g].points[inverse[(*selected)[k]]];continue;
                }
                double best=std::numeric_limits<double>::infinity(); int bestShift=0,bestDirection=1;
                for(int direction:{1,-1}) for(int offset=0;offset<(m_profileClosed?N:1);++offset) {
                    int shift=m_profileClosed?offset:(direction==1?0:N-1);
                    double cost=0;for(int k=0;k<N;++k) cost+=(vertices[stations[g]*N+k]-guides[g].points[(shift+direction*k+N)%N]).squaredNorm();
                    if(cost<best){best=cost;bestShift=shift;bestDirection=direction;}
                }
                for(int k=0;k<N;++k) targets[g][k]=guides[g].points[(bestShift+bestDirection*k+N)%N];
            }
            const int A=std::max(0,M-3), B=std::max(0,M-2);
            residual.resize(3*G*N+3*G+A+3*B+P); double sum=0;maximum=0;int cursor=0;
            double sectionEnergy=0,planeEnergy=0,apexEnergy=0,normalEnergy=0;
            Eigen::MatrixXd J; if(gradient) J=Eigen::MatrixXd::Zero(residual.size(),P);
            std::vector<double> sectionErrors(G,0);
            if(detail) {detail->sectionRms.resize(G);detail->planeAngleDegrees.resize(G);}
            for(int g=0;g<G;++g) {
                int source=guides[g].source;
                double weight=m_settings.targetSectionWeight*(m_settings.sectionWeights.empty()?1:m_settings.sectionWeights[source]);
                double error=0;
                for(int k=0;k<N;++k) {
                    V d=vertices[stations[g]*N+k]-targets[g][k]; error+=d.squaredNorm();maximum=std::max(maximum,d.norm());
                    for(int axis=0;axis<3;++axis) {
                        double factor=std::sqrt(weight/double(N)); double r=factor*d[axis];
                        if(gradient) J.row(cursor)=factor*dv[stations[g]*N+k].row(axis);
                        residual[cursor++]=r;sectionErrors[g]+=r*r;
                    }
                }
                sum+=error;
                if(detail) detail->sectionRms[source]=std::sqrt(error/N);
            }
            // Plane position is fixed by its directrix station; fit the unoriented plane normal.
            for(int g=0;g<G;++g) {
                int source=guides[g].source;
                double weight=m_settings.targetPlaneWeight*(m_settings.planeWeights.empty()?1:m_settings.planeWeights[source]);
                V d=normals[stations[g]].cross(guides[g].normal)*std::sqrt(weight/double(G));
                for(int axis=0;axis<3;++axis) {
                    if(gradient) for(int k=0;k<P;++k) J(cursor,k)=V(dn[stations[g]].col(k)).cross(guides[g].normal)[axis]*std::sqrt(weight/double(G));
                    residual[cursor++]=d[axis];planeEnergy+=d[axis]*d[axis];
                }
                if(detail) detail->planeAngleDegrees[source]=std::acos(std::clamp(std::abs(normals[stations[g]].dot(guides[g].normal)),0.0,1.0))*180.0/3.141592653589793;
            }
            // Mean squared second differences, as in the reference smoothness energies.
            // Apply to actual apex scales, not log controls. Short sequences contribute zero.
            for(int j=1;j<M-2;++j) {
                double d=m_settings.apexScale*(std::exp(q[j-1])-2*std::exp(q[j])+std::exp(q[j+1]));
                double r=d*std::sqrt(m_settings.apexSmoothnessWeight/(2*double(A)));
                if(gradient) {double f=m_settings.apexScale*std::sqrt(m_settings.apexSmoothnessWeight/(2*double(A)));J(cursor,j-1)=f*std::exp(q[j-1]);J(cursor,j)=-2*f*std::exp(q[j]);J(cursor,j+1)=f*std::exp(q[j+1]);}
                residual[cursor++]=r;apexEnergy+=r*r;
            }
            for(int j=1;j<M-1;++j) {
                V d=(normals[j-1]-2*normals[j]+normals[j+1])*std::sqrt(m_settings.planeSmoothnessWeight/(2*double(B)));
                for(int axis=0;axis<3;++axis) {
                    if(gradient) J.row(cursor)=(dn[j-1].row(axis)-2*dn[j].row(axis)+dn[j+1].row(axis))*std::sqrt(m_settings.planeSmoothnessWeight/(2*double(B)));
                    residual[cursor++]=d[axis];normalEnergy+=d[axis]*d[axis];
                }
            }
            for(int k=0;k<P;++k) {if(gradient) J(cursor,k)=1e-5;residual[cursor++]=1e-5*q[k];}
            // Reference smooth_max is the p=10 norm of per-section mean squared errors.
            double largest=*std::max_element(sectionErrors.begin(),sectionErrors.end()), sumP=0;
            if(largest>0) {for(double e:sectionErrors)sumP+=std::pow(e/largest,10);sectionEnergy=largest*std::pow(sumP,.1);}
            if(gradient) {
                *gradient=2*J.bottomRows(residual.size()-3*G*N).transpose()*residual.tail(residual.size()-3*G*N);
                if(sectionEnergy>0) for(int g=0;g<G;++g) *gradient+=2*std::pow(sectionErrors[g]/sectionEnergy,9)*J.middleRows(g*3*N,3*N).transpose()*residual.segment(g*3*N,3*N);
            }
            rms=std::sqrt(sum/(G*N));
            if(detail) {detail->sectionEnergy=sectionEnergy;detail->planeEnergy=planeEnergy;detail->apexSmoothnessEnergy=apexEnergy;detail->planeSmoothnessEnergy=normalEnergy;detail->objective=sectionEnergy+planeEnergy+apexEnergy+normalEnergy+1e-10*q.squaredNorm();}
            return residual.allFinite();
        } catch(const std::exception&) { return false; }
    };
    Eigen::VectorXd residual; double rms,maximum;
    zCTubeFitReport report; report.active=true;
    if(!evaluate(controls,residual,rms,maximum,&report)) throw std::runtime_error("Initial target fit is singular.");
    report.initialRms=rms;report.initialObjective=report.objective;
    auto objective=[&](const Eigen::VectorXd& q,Eigen::VectorXd& gradient) {
        Eigen::VectorXd r;double a,b;zCTubeFitReport d;
        if(!evaluate(q,r,a,b,&d,&gradient) || !gradient.allFinite()) {gradient=Eigen::VectorXd::Zero(P);return std::numeric_limits<double>::infinity();}
        return d.objective;
    };
    // Diagnostic only: production optimisation never uses finite differences.
    if(m_settings.checkFitGradient) {
        Eigen::VectorXd probe=controls;for(int k=0;k<P;++k)probe[k]+=.003*std::sin(double(k+1));
        Eigen::VectorXd analytic;objective(probe,analytic);double worst=0;
        for(int k=0;k<P;++k) {
            const double h=1e-6;auto plus=probe,minus=probe;plus[k]+=h;minus[k]-=h;Eigen::VectorXd unused;
            double numerical=(objective(plus,unused)-objective(minus,unused))/(2*h);
            if(!std::isfinite(numerical) || !analytic.allFinite()) throw std::runtime_error("Gradient check crossed a singular construction.");
            worst=std::max(worst,std::abs(numerical-analytic[k])/std::max({1.0,std::abs(numerical),std::abs(analytic[k])}));
        }
        report.gradientCheckError=worst;
    }
    LBFGSpp::LBFGSBParam<double> solverSettings;
    solverSettings.m=10;solverSettings.epsilon=1e-5;solverSettings.epsilon_rel=0;
    solverSettings.past=1;solverSettings.delta=1e-10;solverSettings.max_iterations=m_settings.fitIterations;solverSettings.max_linesearch=40;
    LBFGSpp::LBFGSBSolver<double> solver(solverSettings);
    Eigen::VectorXd lower=Eigen::VectorXd::Constant(P,-3),upper=Eigen::VectorXd::Constant(P,3);
    double value=report.initialObjective;Eigen::VectorXd best=controls;double bestValue=value;
    auto trackedObjective=[&](const Eigen::VectorXd& q,Eigen::VectorXd& gradient) {
        ++report.evaluations;double f=objective(q,gradient);if(f<bestValue){bestValue=f;best=q;}return f;
    };
    try {report.iterations=solver.minimize(trackedObjective,controls,value,lower,upper);report.termination=report.iterations>=m_settings.fitIterations?"Iteration limit":"Solver tolerance";}
    catch(const std::runtime_error& e) {report.termination=e.what();}
    controls=best;
    Eigen::VectorXd finalGradient;objective(controls,finalGradient);
    report.gradientNorm=0;for(int k=0;k<P;++k)report.gradientNorm=std::max(report.gradientNorm,std::abs(controls[k]-std::clamp(controls[k]-finalGradient[k],lower[k],upper[k])));
    if(report.termination=="Solver tolerance") report.termination=report.gradientNorm<=solverSettings.epsilon?"Projected gradient tolerance":"Objective change tolerance";
    std::vector<V> vertices,normals; sweep(controls,vertices,&normals);
    zPointArray points;for(auto p:vertices) points.push_back(point(p));
    sectionMesh(m_tube,points,m_profileEdges,M,m_profileBranched);
    zPointArray curve;zIntArray edges;
    for(int j=0;j<M;++j) {
        curve.push_back(point(X[j]));if(j)edges.insert(edges.end(),{j-1,j});
        V x=unit(vertices[j*N+1]-vertices[j*N]), normal=unit(normals[j]), y=unit(normal.cross(x));
        m_frames.push_back({point(X[j]),point(x),point(y),point(normal)});
    }
    zFnGraph(m_sampledDirectrix).create(curve,edges);
    evaluate(controls,residual,rms,maximum,&report);
    report.rms=rms;report.maxError=maximum;report.converged=maximum<=m_settings.fitTolerance;
    report.stations.resize(G); for(int g=0;g<G;++g) report.stations[guides[g].source]=stations[g];
    for(int j=0;j<M-1;++j)report.apexScales.push_back(m_settings.apexScale*std::exp(controls[j]));
    m_fit=std::move(report);m_vertices=std::move(points);m_stations=M;
}

void zTsCTubes::computeStrips()
{
    if(m_vertices.empty()) throw std::logic_error("Compute the tube first.");
    zObjectMeshArray strips; int N=int(m_profile.size());
    const int E=int(m_profileEdges.size()/2);
    std::vector<std::vector<int>> incident(N);
    std::vector<int> parent(E);for(int e=0;e<E;++e){parent[e]=e;incident[m_profileEdges[2*e]].push_back(e);incident[m_profileEdges[2*e+1]].push_back(e);}
    std::function<int(int)> root=[&](int e){return parent[e]==e?e:parent[e]=root(parent[e]);};
    auto other=[&](int e,int v){return m_profileEdges[2*e]==v?m_profileEdges[2*e+1]:m_profileEdges[2*e];};
    // Join opposite collinear edges, including across a branched junction.
    // Verify the entire sweep so a real crease at any station remains a seam.
    for(int v=0;v<N;++v){
        struct Pair{int a,b;double error;};std::vector<Pair> candidates;
        for(size_t a=0;a<incident[v].size();++a)for(size_t b=a+1;b<incident[v].size();++b){
            int ea=incident[v][a],eb=incident[v][b];bool straight=true;double worst=0;
            for(int j=0;j<m_stations;++j){V center=vec(m_vertices[j*N+v]),pa=vec(m_vertices[j*N+other(ea,v)]),pb=vec(m_vertices[j*N+other(eb,v)]);
                V da=pa-center,db=pb-center;double la=da.norm(),lb=db.norm();
                double tolerance=1e-3*std::min(la,lb)+4*std::numeric_limits<float>::epsilon()*std::max({center.cwiseAbs().maxCoeff(),pa.cwiseAbs().maxCoeff(),pb.cwiseAbs().maxCoeff()});
                double distance=da.cross(db).norm()/(la+lb);
                if(la<1e-10||lb<1e-10||da.dot(db)>=0||distance>tolerance){straight=false;break;}
                worst=std::max(worst,distance/std::max(tolerance,1e-12));
            }
            if(straight)candidates.push_back({ea,eb,worst});
        }
        std::sort(candidates.begin(),candidates.end(),[](const Pair& a,const Pair& b){return a.error<b.error;});
        std::set<int> paired;for(auto c:candidates)if(!paired.count(c.a)&&!paired.count(c.b)){parent[root(c.b)]=root(c.a);paired.insert(c.a);paired.insert(c.b);}
    }
    std::vector<std::vector<int>> groups(E);for(int e=0;e<E;++e)groups[root(e)].push_back(e);
    double stripExtent=0;for(const auto& p:m_vertices)stripExtent=std::max(stripExtent,(vec(p)-vec(m_vertices.front())).norm());
    for(const auto& group:groups)if(!group.empty())
    {
        std::vector<std::vector<int>> neighbors(N);
        for(int e:group){int a=m_profileEdges[2*e],b=m_profileEdges[2*e+1];neighbors[a].push_back(b);neighbors[b].push_back(a);}
        int start=m_profileEdges[2*group.front()];for(int v=0;v<N;++v)if(neighbors[v].size()==1){start=v;break;}
        zIntArray path;int current=start,previous=-1;
        while(current>=0){path.push_back(current);int next=-1;for(int v:neighbors[current])if(v!=previous){next=v;break;}previous=current;current=next;}
        // A straight flange is one face column. Retaining its internal T anchor
        // creates redundant unfolding cycles, especially at narrow stations.
        // Remove only vertices lying between their neighbours at every section.
        for(size_t k=1;k+1<path.size();) {
            bool redundant=true;
            for(int j=0;j<m_stations;++j) {
                V a=vec(m_vertices[j*N+path[k-1]]),b=vec(m_vertices[j*N+path[k]]),c=vec(m_vertices[j*N+path[k+1]]),ac=c-a;
                double t=(b-a).dot(ac)/ac.squaredNorm();
                double eps=8*std::numeric_limits<float>::epsilon()*std::max({1.0,a.norm(),b.norm(),c.norm()});
                if(!std::isfinite(t)||t<=0||t>=1||(b-a-t*ac).norm()>eps){redundant=false;break;}
                if(j) {
                    V pa=vec(m_vertices[(j-1)*N+path[k-1]]),pc=vec(m_vertices[(j-1)*N+path[k+1]]);
                    V normal=(a-pa).cross(c-pa);
                    // Match the kernel's default planar-face precision; thin
                    // panels magnify rounding of the float vertex storage.
                    if(normal.norm()<eps*eps||std::abs((pc-pa).dot(normal.normalized()))>1e-5*std::max(stripExtent,1e-6)){redundant=false;break;}
                }
            }
            if(redundant)path.erase(path.begin()+k);else ++k;
        }
        const int width=int(path.size());
        zPointArray p; zIntArray counts,connects;
        for(int j=0;j<m_stations;++j)
        {
            for(int v:path)p.push_back(m_vertices[j*N+v]);
            if(j)for(int k=0;k<width-1;++k){int a=width*(j-1)+k;counts.push_back(4);connects.insert(connects.end(),{a,a+width,a+width+1,a+1});}
        }
        strips.emplace_back(); mesh(strips.back(),p,counts,connects);
    }
    m_strips=std::move(strips); m_unrolled.clear(); m_stripAreas.clear();
}
void zTsCTubes::unrollStrips(double gap)
{
    if(!std::isfinite(gap) || gap<0) throw std::invalid_argument("Strip gap must be nonnegative.");
    if(m_strips.empty()) computeStrips();
    zObjectMeshArray flat; double offset=0;
    std::vector<zCTubeStripArea> areas;
    auto area=[](zObjectMesh input) {
        zFnMesh fn(input); fn.triangulate(); zDoubleArray values;
        return fn.getPlanarFaceAreas(values);
    };
    for(auto& strip:m_strips)
    {
        flat.emplace_back();
        try { zFnMesh(strip).unroll(flat.back()); }
        catch(const std::invalid_argument& error) {
            throw std::invalid_argument("Strip "+std::to_string(flat.size())+" cannot unroll exactly: "+error.what());
        }
        zFnMesh fn(flat.back()); zPointArray p; fn.getVertexPositions(p);
        float minX=p[0].x,minY=p[0].y,maxY=p[0].y;
        for(auto& v:p) { minX=std::min(minX,v.x); minY=std::min(minY,v.y); maxY=std::max(maxY,v.y); }
        for(auto& v:p) { v.x-=minX; v.y=float(v.y-minY+offset); }
        fn.setVertexPositions(p); offset+=maxY-minY+gap;
        double a=area(strip), b=area(flat.back());
        if(!std::isfinite(a) || !std::isfinite(b) || a<=0) throw std::runtime_error("Invalid strip area.");
        areas.push_back({a,b,b-a,100.0*std::abs(b-a)/a});
    }
    m_unrolled=std::move(flat); m_stripAreas=std::move(areas);
}
}
