#include <zspace/zToolsets/statics/zTsForceCTubes.h>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace zSpace {
namespace {
using V=Eigen::Vector3d;
using P=Eigen::Vector2d;
V v(const zPoint& p){return {p.x,p.y,p.z};}
zPoint point(const V& p){return {float(p.x()),float(p.y()),float(p.z())};}
void require(bool value,const char* message){if(!value)throw std::invalid_argument(message);}
void finite(const zPointArray& points){for(auto p:points)require(v(p).allFinite()&&v(p).cwiseAbs().maxCoeff()<1e12,"Invalid tube coordinates.");}
struct Properties {double area=0;P c=P::Zero();Eigen::Matrix2d covariance=Eigen::Matrix2d::Zero();std::vector<P> boundary;};
double cross(const P& a,const P& b){return a.x()*b.y()-a.y()*b.x();}
Properties properties(const std::vector<P>& p,const zIntArray& edges,zForceCTubeMaterial mode,double thickness){
    Properties out;Eigen::Matrix2d raw=Eigen::Matrix2d::Zero();P first=P::Zero();
    if(mode==zForceCTubeMaterial::ThinWall){
        for(size_t i=0;i<edges.size();i+=2){
            P a=p[edges[i]],b=p[edges[i+1]],d=b-a,c=(a+b)/2;double length=d.norm();
            require(length>1e-10,"Collapsed section edge.");P u=d/length,n(-u.y(),u.x());double area=length*thickness;
            out.area+=area;first+=area*c;raw+=area*(c*c.transpose()+length*length/12*u*u.transpose()+thickness*thickness/12*n*n.transpose());
            out.boundary.push_back(a+n*thickness/2);out.boundary.push_back(a-n*thickness/2);
            out.boundary.push_back(b+n*thickness/2);out.boundary.push_back(b-n*thickness/2);
        }
    }else{
        // Trace topology: the vertex buffer need not be in polygon order.
        std::vector<std::vector<int>> adj(p.size());
        for(size_t i=0;i<edges.size();i+=2){adj[edges[i]].push_back(edges[i+1]);adj[edges[i+1]].push_back(edges[i]);}
        for(auto& a:adj)require(a.size()==2,"Solid material requires one closed polygon; use Thin wall for an open profile.");
        std::vector<int> order;int previous=-1,current=0;
        do {require(order.size()<p.size(),"Disconnected material polygon.");order.push_back(current);int next=adj[current][0]==previous?adj[current][1]:adj[current][0];previous=current;current=next;}while(current!=0);
        require(order.size()==p.size(),"Disconnected material polygon.");
        const int count=int(order.size());
        for(int i=0;i<count;++i)for(int j=i+1;j<count;++j){
            if(j==i+1||(i==0&&j==count-1))continue;
            P a=p[order[i]],b=p[order[(i+1)%count]],c=p[order[j]],d=p[order[(j+1)%count]];
            double ab1=cross(b-a,c-a),ab2=cross(b-a,d-a),cd1=cross(d-c,a-c),cd2=cross(d-c,b-c);
            require(!(ab1*ab2<0&&cd1*cd2<0),"Self-intersecting material polygon.");
        }
        for(int i=0;i<count;++i){
            P a=p[order[i]],b=p[order[(i+1)%count]];double k=cross(a,b);
            out.area+=k/2;first+=k*(a+b)/6;
            raw(0,0)+=k*(a.x()*a.x()+a.x()*b.x()+b.x()*b.x())/12;
            raw(1,1)+=k*(a.y()*a.y()+a.y()*b.y()+b.y()*b.y())/12;
            raw(0,1)+=k*(2*a.x()*a.y()+a.x()*b.y()+b.x()*a.y()+2*b.x()*b.y())/24;
            out.boundary.push_back(a);
        }
        raw(1,0)=raw(0,1);if(out.area<0){out.area=-out.area;first=-first;raw=-raw;}
    }
    require(std::isfinite(out.area)&&out.area>1e-12,"Section material area is zero or invalid.");
    out.c=first/out.area;out.covariance=raw-out.area*out.c*out.c.transpose();
    require(out.covariance.allFinite()&&out.covariance.determinant()>1e-20,"Degenerate section inertia.");return out;
}
}
struct zTsForceCTubes::Impl {
    zObjectGraph guide,profile,form;zPointArray guidePoints,profilePoints,formPoints;
    zIntArray edges,profileEdges;zDoubleArray areas,forces,differences,apex;
    zForceCTubeSettings settings;zForceCTubeReport report;
    zTsCTubes tube;std::vector<zCTubeFrame> frames;std::vector<zForceCTubeSection> sections;
    bool ready=false,hasForces=false,sourceEquilibrium=false,branched=false;
};
zTsForceCTubes::zTsForceCTubes():impl(new Impl){}
zTsForceCTubes::~zTsForceCTubes()=default;
void zTsForceCTubes::setInput(zObjectGraph& directrix,zObjectGraph& generatrix){
    Impl next;next.settings=impl->settings;next.guide=directrix;next.profile=generatrix;
    zFnGraph(next.guide).getVertexPositions(next.guidePoints);zFnGraph(next.guide).getEdgeData(next.edges);
    zFnGraph(next.profile).getVertexPositions(next.profilePoints);zFnGraph(next.profile).getEdgeData(next.profileEdges);
    const int n=int(next.guidePoints.size()),m=int(next.profilePoints.size());
    require(n>=2&&n<=128&&m>=3&&m<=128,"Force C-Tubes supports 2..128 guide stations and 3..128 profile vertices.");
    require(next.edges.size()==size_t(2*(n-1)),"Force C-Tubes requires one open chain.");
    std::vector<bool> seen(n-1,false);
    for(size_t i=0;i<next.edges.size();i+=2){int a=next.edges[i],b=next.edges[i+1];require(a>=0&&b>=0&&a<n&&b<n&&std::abs(a-b)==1,"Guide vertices must follow chain order.");int e=std::min(a,b);require(!seen[e],"Duplicate guide edge.");seen[e]=true;}
    finite(next.guidePoints);finite(next.profilePoints);next.ready=true;
    // zTsCTubes orders simple sections by traversal and splits branched strip seams.
    // Normalise simple input here so material properties use the same vertex IDs.
    std::vector<std::vector<int>> adj(m);
    require(next.profileEdges.size()%2==0,"Invalid profile edges.");
    for(size_t e=0;e<next.profileEdges.size();e+=2){int a=next.profileEdges[e],b=next.profileEdges[e+1];require(a>=0&&b>=0&&a<m&&b<m,"Invalid profile vertex.");adj[a].push_back(b);adj[b].push_back(a);}
    for(const auto& a:adj){require(!a.empty(),"Disconnected profile vertex.");next.branched=next.branched||a.size()>2;}
    if(!next.branched){int start=0;for(int i=0;i<m;++i)if(adj[i].size()==1){start=i;break;}zPointArray ordered;zIntArray edges;std::vector<bool> visited(m,false);int previous=-1,current=start;
        while(!visited[current]){visited[current]=true;ordered.push_back(next.profilePoints[current]);int candidate=-1;for(int other:adj[current])if(other!=previous){candidate=other;break;}previous=current;if(candidate<0)break;current=candidate;}
        require(ordered.size()==size_t(m),"Disconnected profile.");for(int i=1;i<m;++i)edges.insert(edges.end(),{i-1,i});if(adj[start].size()==2)edges.insert(edges.end(),{m-1,0});next.profilePoints=ordered;next.profileEdges=edges;zFnGraph(next.profile).create(ordered,edges);
    }
    impl.reset(new Impl(std::move(next)));
}
void zTsForceCTubes::setSettings(const zForceCTubeSettings& s){
    for(double x:{s.thickness,s.forceScale,s.referenceStress,s.areaTolerance,s.correspondingAreaTolerance,s.smoothness,s.minScale,s.maxScale})require(std::isfinite(x),"Tube settings must be finite.");
    require(s.thickness>0&&s.forceScale>0&&s.referenceStress>0&&s.areaTolerance>=0&&s.correspondingAreaTolerance>=0&&s.smoothness>=0&&s.minScale>0&&s.maxScale>=s.minScale&&s.minScale<=1&&s.maxScale>=1&&s.maxIterations>=0&&s.maxIterations<=200,"Invalid Force C-Tube settings.");impl->settings=s;
}
void zTsForceCTubes::setForceData(zObjectGraph& form,const zDoubleArray& areas,const zDoubleArray& forces,const zDoubleArray& differences,bool equilibrium){
    auto& t=*impl;require(t.ready,"Set tube inputs first.");zPointArray p;zIntArray edges;zFnGraph(form).getVertexPositions(p);zFnGraph(form).getEdgeData(edges);
    require(p.size()==t.guidePoints.size()&&edges==t.edges,"Force and guide graphs must have matching indexed topology.");
    require(areas.size()*2==edges.size()&&forces.size()==areas.size()&&differences.size()==areas.size(),"Supply one force, face area and paired difference per guide edge.");finite(p);
    for(size_t i=0;i<areas.size();++i)require(std::isfinite(areas[i])&&areas[i]>1e-12&&std::isfinite(forces[i])&&std::isfinite(differences[i])&&differences[i]>=0,"Member face areas must be positive; force data must be finite.");
    t.form=form;t.formPoints=p;t.areas=areas;t.forces=forces;t.differences=differences;t.sourceEquilibrium=equilibrium;t.hasForces=true;
}
void zTsForceCTubes::setForceDiagram(const zTsForcePolyhedra& force,double scale){
    auto graph=force.resultGraph();zDoubleArray sums(force.axialForces().size(),0),counts(sums.size(),0),lo(sums.size(),std::numeric_limits<double>::infinity()),hi(sums.size(),0);
    for(size_t node=0;node<impl->guidePoints.size();++node)for(const auto& face:force.cellFaces(int(node))){
        double total=0;for(auto source:face.sources)total+=source.magnitude;
        for(auto source:face.sources)if(source.type==zForcePolyhedraFaceType::Member&&total>0){int e=source.id;require(e>=0&&e<int(sums.size()),"Invalid force-face source.");double area=face.area*source.magnitude/total;sums[e]+=area;counts[e]++;lo[e]=std::min(lo[e],area);hi[e]=std::max(hi[e],area);}
    }
    zDoubleArray differences(sums.size());for(size_t e=0;e<sums.size();++e){require(counts[e]>0,"Missing member force face.");sums[e]/=counts[e];differences[e]=hi[e]-lo[e];}
    auto settings=impl->settings;settings.forceScale=scale;setSettings(settings);setForceData(graph,sums,force.axialForces(),differences,force.report().equilibrium);
}
zForceCTubeReport zTsForceCTubes::compute(bool optimise){
    auto& t=*impl;require(t.ready&&t.hasForces,"Set tube geometry and a force snapshot first.");const auto& s=t.settings;
    const int M=int(t.guidePoints.size()),N=int(t.profilePoints.size());t.report={};t.report.sourceEquilibrium=t.sourceEquilibrium;
    t.report.forceAreasResolved=std::all_of(t.differences.begin(),t.differences.end(),[&](double d){return d<=s.correspondingAreaTolerance;});
    t.tube.setDirectrix(t.guide);t.tube.setGeneratrix(t.profile);
    // A Rhino section already carries its section plane and roll. Honour a profile
    // vertex on the guide start as the geometric anchor rather than recentering it.
    for(auto p:t.profilePoints)if((v(p)-v(t.guidePoints.front())).norm()<1e-5){t.tube.setGeneratrixOrigin(0,p);break;}
    zCTubeSettings cs;cs.samples=M;cs.spline=false;cs.preserveOrientation=true;
    t.tube.setSettings(cs);t.tube.computeTube();t.frames=t.tube.frames();
    auto sectionVertices=[&](){zPointArray mesh;zFnMesh(t.tube.tubeMesh()).getVertexPositions(mesh);if(!t.branched)return mesh;
        zPointArray sections(size_t(M)*N);for(size_t e=0;e<t.profileEdges.size();e+=2)for(int i=0;i<M;++i){size_t base=e*M+2*i;sections[i*N+t.profileEdges[e]]=mesh[base];sections[i*N+t.profileEdges[e+1]]=mesh[base+1];}return sections;};
    zPointArray seedMesh=sectionVertices();
    zPointArray seed(seedMesh.begin(),seedMesh.begin()+N);
    zDoubleArray target(M,0);std::vector<bool> tension(M,false);
    for(size_t e=0;e<t.areas.size();++e)for(int k=0;k<2;++k){int i=t.edges[2*e+k];target[i]=std::max(target[i],t.areas[e]/s.forceScale/s.referenceStress);tension[i]=tension[i]||t.forces[e]>0;}
    std::vector<Properties> props;zPointArray vertices;
    auto evaluate=[&](const Eigen::VectorXd& q){
        zPointArray scaled;for(auto p:seed)scaled.push_back(point(v(t.frames[0].origin)+std::exp(q[0])*(v(p)-v(t.frames[0].origin))));
        zObjectGraph profile;zFnGraph(profile).create(scaled,t.profileEdges);t.tube.setGeneratrix(profile);
        zDoubleArray apex;for(int i=1;i<M;++i)apex.push_back(std::exp(q[i]));
        t.tube.computeTubeWithFrames(t.frames,{},-1,apex);vertices=sectionVertices();props.clear();
        Eigen::VectorXd r(M+std::max(0,M-2));
        for(int i=0;i<M;++i){std::vector<P> local;for(int k=0;k<N;++k){V d=v(vertices[i*N+k])-v(t.frames[i].origin);local.emplace_back(d.dot(v(t.frames[i].x)),d.dot(v(t.frames[i].y)));}props.push_back(properties(local,t.profileEdges,s.material,s.thickness));r[i]=std::log(props.back().area/target[i]);}
        for(int i=1;i<M-1;++i)r[M+i-1]=std::sqrt(s.smoothness)*(std::log(props[i-1].area)-2*std::log(props[i].area)+std::log(props[i+1].area));return r;
    };
    auto maxError=[&](){double error=0;for(int i=0;i<M;++i)error=std::max(error,std::abs(props[i].area-target[i])/target[i]);return error;};
    Eigen::VectorXd q=Eigen::VectorXd::Zero(M),r=evaluate(q);t.report.initialAreaError=maxError();double lambda=1e-3;
    t.report.termination=optimise?"iteration_limit":"generated";
    if(optimise)for(int iter=0;iter<s.maxIterations;++iter){
        if(maxError()<=s.areaTolerance){t.report.termination="area_tolerance";break;}
        Eigen::MatrixXd J(r.size(),M);const double epsilon=0.001;
        for(int j=0;j<M;++j){auto probe=q;probe[j]+=epsilon;try{J.col(j)=(evaluate(probe)-r)/epsilon;}catch(const std::exception&){J.col(j).setZero();}}
        Eigen::MatrixXd normal=J.transpose()*J;normal.diagonal().array()+=lambda;
        Eigen::VectorXd step=normal.ldlt().solve(-J.transpose()*r);
        if(!step.allFinite()){t.report.termination="singular_step";break;}
        double largest=step.cwiseAbs().maxCoeff();if(largest>0.3)step*=0.3/largest;
        Eigen::VectorXd trial=(q+step).cwiseMax(std::log(s.minScale)).cwiseMin(std::log(s.maxScale));
        bool accepted=false;try{auto next=evaluate(trial);if(next.squaredNorm()<r.squaredNorm()){q=trial;r=next;lambda=std::max(1e-8,lambda/3);accepted=true;}}catch(const std::exception&){}
        t.report.iterations=iter+1;
        if(!accepted){lambda*=10;evaluate(q);if(lambda>1e8){t.report.termination="stalled_or_bounded";break;}}
    }
    evaluate(q);t.apex.clear();for(int i=1;i<M;++i)t.apex.push_back(std::exp(q[i]));
    t.report.maxAreaError=maxError();t.report.areaResolved=t.report.maxAreaError<=s.areaTolerance;t.report.sectionChecksPassed=true;t.sections.clear();
    if(optimise&&t.report.areaResolved)t.report.termination="area_tolerance";
    for(int i=0;i<M;++i){
        const auto& f=t.frames[i];const auto& prop=props[i];zForceCTubeSection out;out.area=prop.area;out.targetArea=target[i];out.relativeAreaError=std::abs(out.area-out.targetArea)/out.targetArea;
        V center=v(f.origin)+prop.c.x()*v(f.x)+prop.c.y()*v(f.y);out.centroid=point(center);
        V hit=V::Zero(),force=V::Zero();bool found=false,ambiguous=false;const V normal=v(f.normal).normalized();
        for(size_t e=0;e<t.areas.size();++e){int a=t.edges[2*e],b=t.edges[2*e+1];V start=v(t.formPoints[a]),d=v(t.formPoints[b])-start;double denominator=normal.dot(d);if(std::abs(denominator)<1e-10*d.norm())continue;
            double u=normal.dot(v(f.origin)-start)/denominator;if(u< -1e-6||u>1+1e-6)continue;V p=start+std::clamp(u,0.0,1.0)*d;
            if(found&&(p-hit).norm()>1e-5){ambiguous=true;continue;}
            // At a load node take the larger one-sided normal-force demand.
            V candidate=d.normalized()*(t.areas[e]/s.forceScale);if(candidate.dot(normal)<0)candidate=-candidate;
            if(!found||candidate.dot(normal)>force.dot(normal)){hit=p;force=candidate;}found=true;
        }
        out.intersectionValid=found&&!ambiguous;
        if(out.intersectionValid){
            out.thrustPoint=point(hit);V offset=hit-center;out.ex=offset.dot(v(f.x));out.ey=offset.dot(v(f.y));out.eccentricity=std::hypot(out.ex,out.ey);
            double axial=force.dot(normal);P eccentric(out.ex,out.ey),gradient=prop.covariance.ldlt().solve(axial*eccentric);
            out.minStress=std::numeric_limits<double>::infinity();out.maxStress=-out.minStress;
            for(const P& p:prop.boundary){double stress=axial/prop.area+gradient.dot(p-prop.c);out.minStress=std::min(out.minStress,stress);out.maxStress=std::max(out.maxStress,stress);}
            out.kernUtilisation=std::max(0.0,1-out.minStress/(axial/prop.area));out.stressUtilisation=out.maxStress/s.referenceStress;
            out.compressionOnly=!tension[i]&&out.minStress>=-1e-7*std::max(1.0,out.maxStress);out.withinStress=out.maxStress<=s.referenceStress*(1+1e-7);
            out.status=tension[i]?"tension_member":!out.compressionOnly?"outside_kern":!out.withinStress?"stress_exceeded":"section_pass";
        }else out.status=ambiguous?"ambiguous_thrust_intersection":"no_thrust_intersection";
        t.report.sectionChecksPassed=t.report.sectionChecksPassed&&out.intersectionValid&&out.compressionOnly&&out.withinStress;
        t.report.maxEccentricity=std::max(t.report.maxEccentricity,out.eccentricity);t.sections.push_back(out);
    }
    for(int i=1;i<M;++i)for(size_t e=0;e<t.profileEdges.size();e+=2){int a=t.profileEdges[e],b=t.profileEdges[e+1];V p=v(vertices[(i-1)*N+a]),u=v(vertices[(i-1)*N+b])-p,w=v(vertices[i*N+a])-p,n=u.cross(w);double scale=std::max(u.norm(),w.norm());if(n.norm()>1e-12&&scale>0)t.report.maxQuadWarp=std::max(t.report.maxQuadWarp,std::abs(n.normalized().dot(v(vertices[i*N+b])-p))/scale);}
    t.report.status=!t.sourceEquilibrium?"provisional_force_snapshot":!t.report.forceAreasResolved?"force_areas_pending":!t.report.areaResolved?"section_areas_pending":!t.report.sectionChecksPassed?"section_checks_pending":"section_checks_pass";return t.report;
}
zObjectMesh& zTsForceCTubes::tubeMesh(){return impl->tube.tubeMesh();}
const std::vector<zCTubeFrame>& zTsForceCTubes::frames()const{return impl->frames;}
const std::vector<zForceCTubeSection>& zTsForceCTubes::sections()const{return impl->sections;}
const zDoubleArray& zTsForceCTubes::apexScales()const{return impl->apex;}
const zForceCTubeReport& zTsForceCTubes::report()const{return impl->report;}
}
