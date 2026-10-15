#include "zTsForcePolyhedraInternal.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace zSpace::forcePolyhedra {
namespace {
double angle(const V& a,const V& b){return std::atan2(a.cross(b).norm(),a.dot(b));}
bool onArc(const V& p,const V& a,const V& b){return std::abs(angle(a,p)+angle(p,b)-angle(a,b))<1e-8;}
V centroid(const Cell& c,const zIntArray& face){V p=V::Zero();for(int i:face)p+=c.positions[i];return p/double(face.size());}
V areaVector(const Cell& c,const zIntArray& f){
    V out=V::Zero(),center=centroid(c,f);
    for(int i=0;i<int(f.size());++i)out+=(c.positions[f[i]]-center).cross(c.positions[f[(i+1)%f.size()]]-center)*0.5;
    return out;
}
double correctionEnergy(const Cell& c,std::vector<V>* gradient=nullptr){
    if(gradient)gradient->assign(c.positions.size(),V::Zero());
    double scale=0,energy=0,extent=std::max(c.dynamics->scale,1e-12);
    for(const auto& face:c.faces)if(!face.sources.empty())scale=std::max(scale,face.targetArea);
    for(int f=0;f<int(c.faces.size());++f){
        const auto& face=c.faces[f];const auto& polygon=c.polygons[f];V a=areaVector(c,polygon),center=centroid(c,polygon);double area=a.norm();
        V n=area>1e-25?a/area:c.normals[f],derivative=V::Zero();
        if(face.type!=zForcePolyhedraFaceType::Virtual){double denominator=face.targetArea>0?face.targetArea:std::max(scale,1e-25),error=(area-face.targetArea)/denominator;energy+=.5*error*error;derivative+=n*(error/denominator);}
        if(face.type==zForcePolyhedraFaceType::Zero){
            Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero();for(int id:polygon){V d=c.positions[id]-center;covariance+=d*d.transpose();}
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);V axis=eigen.eigenvectors().col(2);
            for(int id:polygon){V d=c.positions[id]-center;d-=axis*axis.dot(d);energy+=5*d.squaredNorm()/(extent*extent);if(gradient)(*gradient)[id]+=10*d/(extent*extent);}
        }
        if(face.type!=zForcePolyhedraFaceType::Zero&&area>1e-25){
            double alignment=std::clamp(n.dot(c.normals[f]),-1.0,1.0);energy+=5*(1-alignment);
            derivative-=5*(c.normals[f]-n*alignment)/area;
            V planeDerivative=V::Zero();
            for(int id:polygon){V d=c.positions[id]-center;double distance=d.dot(n);energy+=.5*distance*distance/(extent*extent);planeDerivative+=d*(distance/(extent*extent));if(gradient)(*gradient)[id]+=n*(distance/(extent*extent));}
            V centerDerivative=V::Zero();
            for(int id=0;id<int(c.positions.size());++id){
                V d=c.positions[id]-center;double distance=d.dot(n);if(distance<=0)continue;
                energy+=2.5*distance*distance/(extent*extent);double coefficient=5*distance/(extent*extent);
                planeDerivative+=coefficient*d;centerDerivative-=coefficient*n;if(gradient)(*gradient)[id]+=coefficient*n;
            }
            if(gradient)for(int id:polygon)(*gradient)[id]+=centerDerivative/double(polygon.size());
            derivative+=(planeDerivative-n*n.dot(planeDerivative))/area;
        }
        if(gradient)for(int i=0;i<int(polygon.size());++i){V edge=c.positions[polygon[(i+1)%polygon.size()]]-c.positions[polygon[(i+polygon.size()-1)%polygon.size()]];(*gradient)[polygon[i]]+=.5*edge.cross(derivative);}
    }
    return energy;
}
double acceptIncrement(Cell& c,const std::vector<V>& previous,const std::vector<V>& delta,const zForcePolyhedraSettings& s,double before){
    std::vector<V> normals;for(const auto& polygon:c.polygons)normals.push_back(areaVector(c,polygon));
    const double cosine=std::cos(s.normalStepDegrees*3.14159265358979323846/180.0);
    double alpha=1;
    for(int attempt=0;attempt<22;++attempt,alpha*=.5){
        for(int i=0;i<int(previous.size());++i)c.positions[i]=previous[i]+alpha*delta[i];
        bool valid=true;
        for(int f=0;f<int(normals.size())&&valid;++f)if(c.faces[f].type!=zForcePolyhedraFaceType::Zero&&normals[f].norm()>1e-15){V now=areaVector(c,c.polygons[f]);double denominator=now.norm()*normals[f].norm();valid=denominator>1e-30&&now.dot(normals[f])>=cosine*denominator;}
        if(valid&&correctionEnergy(c)<=before+1e-12*std::max(1e-12,before))return alpha;
    }
    c.positions=previous;return 0;
}
void publish(Cell& c){
    zPointArray positions;zIntArray counts,connects;
    for(const V& p:c.positions){require(p.allFinite()&&p.cwiseAbs().maxCoeff()<1e20,"Nonfinite area-pursuit geometry.");positions.push_back(publicVector(p));}
    for(const auto& f:c.polygons){counts.push_back(int(f.size()));connects.insert(connects.end(),f.begin(),f.end());}
    if(!positions.empty())zFnMesh(c.mesh).create(positions,counts,connects);
}
void measure(Cell& c,const zForcePolyhedraSettings& s){
    auto& r=c.report;r.maxAreaError=r.maxNormalError=r.maxPlanarityError=0;
    double scale=0,extent=0;for(const auto& f:c.faces)if(f.type!=zForcePolyhedraFaceType::Virtual)scale=std::max(scale,f.targetArea);
    for(const V& p:c.positions)extent=std::max(extent,p.norm());
    V closure=V::Zero();for(int i=0;i<int(c.faces.size());++i)if(!c.faces[i].sources.empty())closure+=c.normals[i]*c.faces[i].targetArea;
    r.closureError=closure.norm()/s.forceScale;
    for(int i=0;i<int(c.faces.size());++i){
        V a=areaVector(c,c.polygons[i]),center=centroid(c,c.polygons[i]);auto& f=c.faces[i];f.area=a.norm();
        if(f.type!=zForcePolyhedraFaceType::Virtual)r.maxAreaError=std::max(r.maxAreaError,std::abs(f.area-f.targetArea)/(f.targetArea>0?f.targetArea:scale));
        // Vanishing zero faces have no measurable normal and do not limit convergence.
        if(f.type!=zForcePolyhedraFaceType::Zero&&f.area>scale*s.areaTolerance){
            double alignment=a.dot(c.normals[i])/f.area;
            r.maxNormalError=std::max(r.maxNormalError,alignment<=0?1.0:a.cross(c.normals[i]).norm()/f.area);
            for(int v:c.polygons[i])r.maxPlanarityError=std::max(r.maxPlanarityError,std::abs((c.positions[v]-center).dot(c.normals[i]))/std::max(1e-15,extent));
            // Convexity: every vertex must lie behind every nonvanishing face.
            // Include supporting-plane violations in the plane error so a
            // tangled or inward cell cannot claim convergence from areas alone.
            for(const V& p:c.positions)r.maxPlanarityError=std::max(r.maxPlanarityError,(p-center).dot(c.normals[i])/std::max(1e-15,extent));
        }
    }
    r.converged=r.maxAreaError<=s.areaTolerance&&r.maxNormalError<=s.normalTolerance&&r.maxPlanarityError<=s.normalTolerance;
    r.status=r.converged?"converged":"area_pursuit_incomplete";
}
}
V faceNormal(const Cell& c,int face){
    V n=areaVector(c,c.polygons[face]);return n.norm()>1e-20?n.normalized():c.normals[face];
}
void refresh(Cell& c,const zForcePolyhedraSettings& s){if(!c.positions.empty())measure(c,s);}
// Complete low-valence physical force sets in their equilibrium dimension.
// Virtual caps must not carry the out-of-plane imbalance of a spatial seed.
void alignVirtualCaps(Cell& c, bool closeTriangle){
    std::vector<int> physical,virtuals;int load=-1;
    for(int i=0;i<int(c.faces.size());++i){if(!c.faces[i].sources.empty()){physical.push_back(i);for(auto source:c.faces[i].sources)if(source.type==zForcePolyhedraFaceType::Load)load=i;}else if(c.faces[i].type==zForcePolyhedraFaceType::Virtual)virtuals.push_back(i);}
    if(physical.size()<2||physical.size()>3)return;
    V axis=c.normals[load>=0?load:physical[0]],plane;
    if(physical.size()==2){
        for(int id:physical)if(id!=(load>=0?load:physical[0]))c.normals[id]=-axis;
        plane=axis.unitOrthogonal();
    }else{
        Eigen::Matrix3d n;for(int i=0;i<3;++i)n.row(i)=c.normals[physical[i]].transpose();
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(n,Eigen::ComputeFullV);plane=svd.matrixV().col(2);plane-=axis*plane.dot(axis);
        if(plane.norm()<1e-12)plane=axis.unitOrthogonal();plane.normalize();
        for(int id:physical)if(id!=load){V projected=c.normals[id]-plane*c.normals[id].dot(plane);if(projected.norm()<1e-12)projected=axis.cross(plane);c.normals[id]=projected.normalized();}
        bool reaction=false;for(int id:physical)for(auto source:c.faces[id].sources)reaction|=source.type==zForcePolyhedraFaceType::Reaction;
        if(closeTriangle && !reaction){
            // A three-force node has a closed force triangle. Seed its two free
            // directions from that triangle, preserving the prescribed load
            // direction and all target magnitudes before the EGI/AP stage.
            int anchor=load>=0?load:physical[0];std::vector<int> other;for(int id:physical)if(id!=anchor)other.push_back(id);
            double p=c.faces[anchor].targetArea,f=c.faces[other[0]].targetArea,g=c.faces[other[1]].targetArea;
            if(p>0&&f>0&&g>0){
                double vertical=(p*p+f*f-g*g)/(2*p),horizontal2=f*f-vertical*vertical;
                if(horizontal2>1e-16*std::max(f*f,g*g)){
                    V sideways=c.normals[other[0]]-axis*c.normals[other[0]].dot(axis);if(sideways.norm()<1e-12)sideways=axis.cross(plane);sideways.normalize();
                    V first=-axis*vertical+sideways*std::sqrt(horizontal2),second=-axis*p-first;
                    c.normals[other[0]]=first.normalized();c.normals[other[1]]=second.normalized();
                }
            }
        }
    }
    if(!virtuals.empty()&&plane.dot(c.normals[virtuals[0]])<0)plane=-plane;
    for(int id:physical)c.faces[id].normal=publicVector(c.normals[id]);
    if(virtuals.size()==2){for(int i=0;i<2;++i){c.normals[virtuals[i]]=plane*(i? -1.0:1.0);c.faces[virtuals[i]].normal=publicVector(c.normals[virtuals[i]]);}}
    if(virtuals.size()==4){V other=axis.cross(plane);for(int i=0;i<4;++i){c.normals[virtuals[i]]=(i<2?plane:other)*(i%2?-1.0:1.0);c.faces[virtuals[i]].normal=publicVector(c.normals[virtuals[i]]);}}
}
void buildCell(Cell& c,const zForcePolyhedraSettings& s,bool allowUnbalanced){
    if(c.normals.empty()){c.report.status="no_forces";c.report.converged=true;return;}
    if(allowUnbalanced)alignVirtualCaps(c,true);
    require(c.normals.size()<=16,"A force cell supports at most 16 distinct physical force directions.");
    Eigen::MatrixXd n(c.normals.size(),3);V closure=V::Zero();double scale=0;
    for(int i=0;i<n.rows();++i){n.row(i)=c.normals[i].transpose();closure+=c.normals[i]*c.faces[i].targetArea;scale=std::max(scale,c.faces[i].targetArea);}
    c.report.closureError=closure.norm()/s.forceScale;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(n,Eigen::ComputeFullV);svd.setThreshold(1e-7);c.report.forceRank=int(svd.rank());
    // An unbalanced fit may satisfy the user's absolute equilibrium tolerance but
    // still be too far from a valid EGI at the relative cell scale. Do not hide it.
    if(!allowUnbalanced&&closure.norm()>scale*std::max(s.areaTolerance,1e-7)){
        c.report.status="unbalanced_cell";return;
    }
    auto add=[&](const V& normal,zForcePolyhedraFaceType type){
        for(int i=0;i<int(c.normals.size());++i)if((normal-c.normals[i]).norm()<1e-8)return i;
        require(c.normals.size()<2048,"EGI arrangement exceeds the cell complexity limit.");
        c.normals.push_back(normal);zForcePolyhedraFace face;face.normal=publicVector(normal);face.type=type;c.faces.push_back(face);return int(c.normals.size()-1);
    };
    // Temporarily unbalanced cells: complete the empty hemisphere with a
    // virtual closing normal, targeted to zero area (paper Fig. 6c). AP may
    // rotate member faces while this temporary face collapses.
    if(allowUnbalanced&&closure.norm()>scale*1e-9)add(-closure.normalized(),zForcePolyhedraFaceType::Zero);
    // Fig. 6: virtual caps have no physical area targets. For rank one add four
    // transverse normals; for rank two add the pair perpendicular to the node.
    // Two or three physical forces must become collinear or coplanar at
    // equilibrium, even if the initial spatial graph is not. Supply caps for
    // that limiting cell before adaptive AP collapses its temporary closing
    // face. Otherwise a three-force tetrahedral seed collapses its whole volume.
    const int capRank=std::min(c.report.forceRank,std::max(1,int(n.rows())-1));
    V fixedLoad=c.normals[0];for(int i=0;i<n.rows();++i)for(auto source:c.faces[i].sources)if(source.type==zForcePolyhedraFaceType::Load)fixedLoad=c.normals[i];
    V transverse=fixedLoad.unitOrthogonal();
    for(int k=capRank;k<3;++k){
        V v=capRank==1?(k==1?transverse:fixedLoad.cross(transverse)):svd.matrixV().col(k);
        if(capRank==2){v-=fixedLoad*v.dot(fixedLoad);if(v.norm()<1e-12)v=fixedLoad.unitOrthogonal();v.normalize();}
        add(v,zForcePolyhedraFaceType::Virtual);add(-v,zForcePolyhedraFaceType::Virtual);
    }
    struct Arc{int u,v;};std::vector<Arc> arcs;int count=int(c.normals.size());
    for(int i=0;i<count;++i)for(int j=i+1;j<count;++j)if((c.normals[i]+c.normals[j]).norm()>1e-9)arcs.push_back({i,j});
    // Complete minor-arc arrangement. Cross-adjacencies become zero EGI masses
    // before extracting spherical regions and taking the topological dual.
    for(int i=0;i<int(arcs.size());++i)for(int j=i+1;j<int(arcs.size());++j){
        auto a=arcs[i],b=arcs[j];V planeA=c.normals[a.u].cross(c.normals[a.v]),planeB=c.normals[b.u].cross(c.normals[b.v]);
        if(planeA.norm()<1e-12||planeB.norm()<1e-12)continue;
        V cross=planeA.normalized().cross(planeB.normalized());
        if(cross.norm()<1e-10)continue;cross.normalize();
        for(double sign:{1.0,-1.0}){V p=sign*cross;if(onArc(p,c.normals[a.u],c.normals[a.v])&&onArc(p,c.normals[b.u],c.normals[b.v]))add(p,zForcePolyhedraFaceType::Zero);}
    }
    count=int(c.normals.size());std::vector<zIntArray> neighbors(count);
    std::set<std::pair<int,int>> edges;
    for(auto arc:arcs){
        std::vector<std::pair<double,int>> path;
        for(int i=0;i<count;++i)if(onArc(c.normals[i],c.normals[arc.u],c.normals[arc.v]))path.emplace_back(angle(c.normals[arc.u],c.normals[i]),i);
        std::sort(path.begin(),path.end());
        for(int i=1;i<int(path.size());++i){int u=path[i-1].second,v=path[i].second;if(edges.emplace(std::min(u,v),std::max(u,v)).second){neighbors[u].push_back(v);neighbors[v].push_back(u);}}
    }
    for(int i=0;i<count;++i){
        V x=c.normals[i].unitOrthogonal(),y=c.normals[i].cross(x);
        auto azimuth=[&](int j){return std::atan2(c.normals[j].dot(y),c.normals[j].dot(x));};
        std::sort(neighbors[i].begin(),neighbors[i].end(),[&](int a,int b){return azimuth(a)<azimuth(b);});
        if(neighbors[i].size()<3)throw std::runtime_error("Invalid EGI adjacency at normal "+std::to_string(i)+" of "+std::to_string(count)+" (degree "+std::to_string(neighbors[i].size())+").");
    }
    std::set<std::pair<int,int>> visited;std::vector<zIntArray> regions;
    for(auto edge:edges)for(auto directed:{edge,std::pair<int,int>(edge.second,edge.first)}){
        if(visited.count(directed))continue;int u=directed.first,v=directed.second;zIntArray region;
        do{
            require(visited.emplace(u,v).second,"Invalid EGI face traversal.");region.push_back(u);
            auto& adjacent=neighbors[v];auto at=std::find(adjacent.begin(),adjacent.end(),u);
            int index=int(at-adjacent.begin());int w=adjacent[(index+int(adjacent.size())-1)%adjacent.size()];u=v;v=w;
            require(region.size()<=edges.size()*2,"EGI traversal did not close.");
        }while(std::pair<int,int>(u,v)!=directed);
        require(region.size()>=3,"Invalid EGI spherical region.");regions.push_back(region);
    }
    require(int(regions.size())-int(edges.size())+count==2,"EGI is not a closed spherical mesh.");
    c.polygons.resize(count);double radius=std::sqrt(scale);
    for(int i=0;i<int(regions.size());++i){
        V p=V::Zero();for(int j:regions[i]){p+=c.normals[j];c.polygons[j].push_back(i);}
        require(p.norm()>1e-10,"EGI region has no finite primal seed.");
        Eigen::MatrixXd planes(regions[i].size(),3);for(int j=0;j<int(regions[i].size());++j)planes.row(j)=c.normals[regions[i][j]].transpose();
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(planes);
        V vertex=solver.rank()==3?V(solver.solve(Eigen::VectorXd::Ones(planes.rows()))):p.normalized();
        c.positions.push_back(vertex*radius);
    }
    // Order dual vertices around each target normal, keeping outward winding.
    for(int i=0;i<count;++i){
        auto& polygon=c.polygons[i];V center=centroid(c,polygon),x=c.normals[i].unitOrthogonal(),y=c.normals[i].cross(x);
        auto azimuth=[&](int j){V p=c.positions[j]-center;return std::atan2(p.dot(y),p.dot(x));};
        std::sort(polygon.begin(),polygon.end(),[&](int a,int b){return azimuth(a)<azimuth(b);});
    }
    double current=0,target=0;for(int i=0;i<int(c.faces.size());++i)if(!c.faces[i].sources.empty()){current+=areaVector(c,c.polygons[i]).norm();target+=c.faces[i].targetArea;}
    if(current>1e-25&&target>0){double factor=std::sqrt(target/current);for(V& p:c.positions)p*=factor;}
    V center=V::Zero();for(const V& p:c.positions)center+=p;center/=double(c.positions.size());for(V& p:c.positions)p-=center;
    measure(c,s);publish(c);
}
void pursue(Cell& c,const zForcePolyhedraSettings& s,int iterations,bool adaptiveNormals){
    if(c.positions.empty())return;
    if(!c.dynamics){
        c.dynamics=std::make_unique<MeshDynamicsState>();
        auto& state=*c.dynamics;double extent=0;for(const V& p:c.positions)extent=std::max(extent,p.norm());
        state.scale=std::max(extent,1e-12);zPointArray points;zIntArray counts,connects;
        for(const V& p:c.positions)points.push_back(publicVector(p/state.scale));
        for(const auto& polygon:c.polygons){counts.push_back(int(polygon.size()));connects.insert(connects.end(),polygon.begin(),polygon.end());}
        zFnMesh(state.object).create(points,counts,connects);state.fn.makeDynamic(false);
    }
    c.report.maxVertexStep=0;
    for(int iteration=0;iteration<iterations&&!c.report.converged;++iteration){
        if(adaptiveNormals)for(int f=0;f<int(c.faces.size());++f){
            bool fixed=c.faces[f].sources.empty();for(auto source:c.faces[f].sources)fixed|=source.type==zForcePolyhedraFaceType::Load;
            if(!fixed){c.normals[f]=faceNormal(c,f);c.faces[f].normal=publicVector(c.normals[f]);}
        }
        std::vector<V> sums(c.positions.size(),V::Zero());std::vector<int> count(c.positions.size(),0);
        for(int f=0;f<int(c.faces.size());++f){
            const auto& polygon=c.polygons[f];const auto& face=c.faces[f];V center=centroid(c,polygon);
            V area=areaVector(c,polygon);double magnitude=area.norm();
            if(face.type==zForcePolyhedraFaceType::Zero){
                Eigen::Matrix3d covariance=Eigen::Matrix3d::Zero();for(int id:polygon){V d=c.positions[id]-center;covariance+=d*d.transpose();}
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance);V axis=eigen.eigenvectors().col(2);
                for(int id:polygon){V d=c.positions[id]-center;sums[id]+=center+axis*axis.dot(d);count[id]++;}
                continue;
            }
            V current=magnitude>1e-25?area/magnitude:c.normals[f],target=c.normals[f];
            V axis=current.cross(target);double rotation=std::atan2(axis.norm(),current.dot(target));
            Eigen::Matrix3d transform=Eigen::Matrix3d::Identity();
            if(rotation>1e-14){if(axis.norm()<1e-12)axis=current.unitOrthogonal();else axis.normalize();transform=Eigen::AngleAxisd(std::min(rotation,s.normalStepDegrees*3.14159265358979323846/180.0),axis).toRotationMatrix();}
            V plane=transform*current;
            double scale=face.type==zForcePolyhedraFaceType::Virtual||magnitude<1e-25?1:std::clamp(std::sqrt(face.targetArea/magnitude),.8,1.25);
            for(int id:polygon){V d=transform*(c.positions[id]-center);d-=plane*plane.dot(d);sums[id]+=center+scale*d;count[id]++;}
        }
        std::vector<V> corrections(c.positions.size());V average=V::Zero();
        for(int i=0;i<int(corrections.size());++i){corrections[i]=sums[i]/count[i]-c.positions[i];average+=corrections[i];}average/=double(corrections.size());
        for(V& correction:corrections)correction-=average;
        auto previous=c.positions;double before=correctionEnergy(c);
        auto delta=integrateCorrections(*c.dynamics,corrections,s,s.maxVertexStepRatio*c.dynamics->scale);++c.report.dynamicsSteps;
        double accepted=acceptIncrement(c,previous,delta,s,before);
        if(!accepted){
            std::vector<V> gradient;correctionEnergy(c,&gradient);
            V mean=V::Zero();for(int i=0;i<int(gradient.size());++i){gradient[i]*=-c.dynamics->scale*c.dynamics->scale/count[i];mean+=gradient[i];}mean/=double(gradient.size());for(V& g:gradient)g-=mean;
            delta=integrateCorrections(*c.dynamics,gradient,s,s.maxVertexStepRatio*c.dynamics->scale);++c.report.dynamicsSteps;
            accepted=acceptIncrement(c,previous,delta,s,before);
        }
        for(const V& step:delta)c.report.maxVertexStep=std::max(c.report.maxVertexStep,accepted*step.norm());
        ++c.report.iterations;measure(c,s);if(!accepted&&!c.report.converged)c.report.status="dynamics_stalled";
    }
    publish(c);
}
}
