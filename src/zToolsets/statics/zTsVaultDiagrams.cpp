#include "zTsVaultDiagrams.h"
#include <Eigen/SparseQR>
#include <map>
#include <set>
#include <numeric>
#include <cmath>
#include <stdexcept>

namespace zSpace {
namespace {
void check(bool b,const char* message){if(!b)throw std::invalid_argument(message);}
struct Union {
    std::vector<int> p;
    explicit Union(int n):p(n){std::iota(p.begin(),p.end(),0);}
    int root(int a){return p[a]==a?a:p[a]=root(p[a]);}
    void join(int a,int b){p[root(a)]=root(b);}
};
using Pair=std::pair<int,int>;
}
zVaultDiagramData buildVaultDiagrams(zObjectMesh& result,const zIntArray& supports,
    const zIntArray& edges,const zDoubleArray& q,const zVaultDiagramSettings& settings) {
    check(std::isfinite(settings.tolerance)&&settings.tolerance>0,"Invalid diagram tolerance.");
    zVaultDiagramData out;out.settings=settings;
    zPointArray positions;zIntArray counts,connects;
    zFnMesh fn(result);fn.getVertexPositions(positions);fn.getPolygonData(connects,counts);
    const int nv=int(positions.size()), ne=int(edges.size()/2), nf=int(counts.size());
    check(q.size()==size_t(ne),"Missing fitted force densities.");
    std::vector<bool> fixed(nv,false);for(int v:supports)fixed[v]=true;
    std::map<Pair,int> face;
    int offset=0;
    for(int f=0;f<nf;++f){
        std::set<int> unique;
        for(int j=0;j<counts[f];++j){
            int a=connects[offset+j],b=connects[offset+(j+1)%counts[f]];
            check(unique.insert(a).second,"Repeated face vertex is unsupported.");
            check(face.emplace(Pair(a,b),f).second,"Diagrams require consistently oriented manifold faces.");
        }
        offset+=counts[f];
    }
    Union exterior(ne);
    std::vector<zIntArray> boundary(nv);
    for(int e=0;e<ne;++e){
        int a=edges[2*e],b=edges[2*e+1];
        if(!face.count({a,b})||!face.count({b,a})){
            boundary[a].push_back(e);boundary[b].push_back(e);
        }
    }
    for(int v=0;v<nv;++v)if(!boundary[v].empty()){
        check(boundary[v].size()==2,"Non-manifold boundary vertex.");
        // A loop without supports (e.g. an opening) is one exterior dual region.
        // Supports split regions; they are not required on every boundary loop.
        if(!fixed[v])exterior.join(boundary[v][0],boundary[v][1]);
    }
    std::map<int,int> exteriorNodes;
    std::map<Pair,int> side=face;
    int nodes=nf;
    for(int e=0;e<ne;++e){
        int a=edges[2*e],b=edges[2*e+1];
        if(fixed[a]&&fixed[b])continue;
        for(Pair pair:{Pair(a,b),Pair(b,a)})if(!side.count(pair)){
            int root=exterior.root(e);
            if(!exteriorNodes.count(root))exteriorNodes[root]=nodes++;
            side[pair]=exteriorNodes[root];
        }
    }
    for(auto& p:positions)p.z=0;
    zFnMesh(out.form).create(positions,counts,connects);
    zIntArray forceCounts,forceConnects;
    // Ordered vertex fans are a topology operation: use zSpace iterators.
    for(zItMeshVertex v(out.form);!v.end();v++)if(!fixed[v.getId()]){
        zItMeshHalfEdgeArray fan;v.getConnectedHalfEdges(fan);
        check(fan.size()>=3,"Stage 1 requires at least three edges at each free vertex.");
        std::set<int> polygon;
        for(auto& h:fan){
            Pair pair(h.getStartVertex().getId(),h.getVertex().getId());
            check(side.count(pair)!=0,"Incomplete dual boundary mapping.");
            int node=side.at(pair);
            check(polygon.insert(node).second,"Unsupported pinched vertex fan in force diagram.");
            forceConnects.push_back(node);
        }
        forceCounts.push_back(int(fan.size()));
    }
    struct Member {int e,a,b;Eigen::Vector2d target;};
    std::vector<Member> members;Union components(nodes);
    Eigen::MatrixXd balance=Eigen::MatrixXd::Zero(nv,2);
    for(int e=0;e<ne;++e){
        int u=edges[2*e],v=edges[2*e+1];if(fixed[u]&&fixed[v])continue;
        Eigen::Vector2d d(positions[v].x-positions[u].x,positions[v].y-positions[u].y);
        check(d.norm()>1e-10,"Vertical or collapsed projected edge is unsupported.");
        int a=side.at({u,v}),b=side.at({v,u});
        check(a!=b,"Force edge has identical dual endpoints.");
        Eigen::Vector2d target=q[e]*d;
        members.push_back({e,a,b,target});components.join(a,b);
        balance.row(u)+=target.transpose();balance.row(v)-=target.transpose();
    }
    std::set<int> anchors;for(int i=0;i<nodes;++i)anchors.insert(components.root(i));
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::MatrixXd rhs=Eigen::MatrixXd::Zero(members.size()+anchors.size(),2);
    int row=0;
    for(const auto& m:members){triplets.emplace_back(row,m.a,-1);triplets.emplace_back(row,m.b,1);rhs.row(row++)=m.target.transpose();}
    for(int a:anchors)triplets.emplace_back(row++,a,1);
    Eigen::SparseMatrix<double> A(row,nodes);A.setFromTriplets(triplets.begin(),triplets.end());
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> qr;
    qr.compute(A);check(qr.info()==Eigen::Success&&qr.rank()==nodes,"Singular reciprocal reconstruction.");
    Eigen::MatrixXd xy=qr.solve(rhs);check(qr.info()==Eigen::Success&&xy.allFinite(),"Reciprocal solve failed.");
    zPointArray forcePositions;for(int i=0;i<nodes;++i)forcePositions.emplace_back(float(xy(i,0)),float(xy(i,1)),0.0f);
    zFnMesh(out.force).create(forcePositions,forceCounts,forceConnects);
    std::map<Pair,int> forceHalfedges;
    zItMeshHalfEdge fh(out.force);out.forceToForm.assign(fh.size(),-1);
    for(;!fh.end();fh++)forceHalfedges[{fh.getStartVertex().getId(),fh.getVertex().getId()}]=fh.getId();
    zItMeshHalfEdge h(out.form);out.formToForce.assign(h.size(),-1);
    for(;!h.end();h++){
        int u=h.getStartVertex().getId(),v=h.getVertex().getId();if(fixed[u]&&fixed[v])continue;
        Pair pair(side.at({u,v}),side.at({v,u}));
        check(forceHalfedges.count(pair)!=0,"Missing reciprocal edge in force mesh.");
        int forceId=forceHalfedges.at(pair);
        check(out.forceToForm[forceId]==-1,"Multiple form edges map to one force edge; unsupported topology.");
        out.formToForce[h.getId()]=forceId;out.forceToForm[forceId]=h.getId();
    }
    for(zItMeshEdge edge(out.form);!edge.end();edge++){
        zIntArray vertices;edge.getVertices(vertices);int source=-1;
        check(fn.edgeExists(vertices[0],vertices[1],source),"Missing source edge.");
        const zColor color=fixed[vertices[0]]&&fixed[vertices[1]]?zColor(0,0,0,1):
            q[source]>0?zColor(240.0f/255,0,140.0f/255,1):zColor(0,40.0f/255,240.0f/255,1);
        edge.setColor(color);
        int mapped=out.formToForce[edge.getHalfEdge(0).getId()];
        if(mapped>=0){zItMeshHalfEdge dual(out.force,mapped);dual.getEdge().setColor(color);}
    }
    double sum=0;
    for(int i=0;i<nv;++i)if(!fixed[i]){double error=balance.row(i).norm();sum+=error*error;out.report.closureError=std::max(out.report.closureError,error);}
    out.report.horizontalResidual=std::sqrt(sum);
    for(const auto& m:members){
        // Evaluate the delivered float geometry, not just the double solve.
        Eigen::Vector2d actual(forcePositions[m.b].x-forcePositions[m.a].x,forcePositions[m.b].y-forcePositions[m.a].y);
        out.report.maxEdgeMismatch=std::max(out.report.maxEdgeMismatch,(actual-m.target).norm());
    }
    out.report.exact=out.report.closureError<=settings.tolerance&&out.report.maxEdgeMismatch<=settings.tolerance;
    out.report.status=out.report.exact?"reciprocal":"approximate";
    return out;
}
}
