#include <zspace/zToolsets/statics/zTsVault.h>
#include "zTsVaultDiagrams.h"
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <limits>
#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>

namespace zSpace {
namespace {
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
void require(bool test, const char* message) { if (!test) throw std::invalid_argument(message); }
bool finite(const zVector& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
}
struct zTsVault::Impl {
    zVaultDiagramData diagrams;
    zObjectMesh mesh;
    zPointArray original, reference;
    zIntArray edges, supports, tension, active, free, row;
    zDoubleArray masses, lower, upper, densities, forces;
    zVector gravity = zVector(0,0,-1);
    zVectorArray loads, residuals, reactions;
    zVaultBestFitSettings settings;
    zVaultSolveReport report;
    Mat X, P;
    Vec q, lo, hi;
    struct State { Mat X; Vec q; zVaultSolveReport report; double h; };
    std::deque<State> history;
    State initial;
    double h = 0.1, minLength = 1e-10;
    double edgeReferenceLength = 1;
    double planReferenceLength = 1;
    std::function<bool(const zPointArray&)> geometryValidator;
    bool initialized = false;
    bool largeSystem() const { return active.size()>1024 || double(active.size())*free.size()*3>8000000; }
    void invalidate() { initialized = false; history.clear(); report = {}; diagrams = {}; }
    void ids(const zIntArray& ids, size_t count) {
        std::vector<bool> seen(count,false);
        for (int id : ids) { require(id>=0 && size_t(id)<count && !seen[id], "Invalid or duplicate component ID."); seen[id]=true; }
    }
    Mat matrix(const Mat& x) const {
        Mat a = Mat::Zero(3*free.size(), active.size());
        for (int j=0;j<int(active.size());++j) {
            int e=active[j], u=edges[2*e], v=edges[2*e+1];
            Eigen::RowVector3d d=x.row(u)-x.row(v);
            require(d.allFinite() && d.norm()>minLength, "Collapsed or nonfinite member.");
            for(int k=0;k<3;++k) {
                if(row[u]>=0) a(row[u]+k*free.size(),j)=d[k];
                if(row[v]>=0) a(row[v]+k*free.size(),j)=-d[k];
            }
        }
        return a;
    }
    Vec loadVector() const {
        Vec p(3*free.size());
        for(int i=0;i<int(free.size());++i) for(int k=0;k<3;++k) p[i+k*free.size()]=P(free[i],k);
        return p;
    }
    bool validGeometry(const Mat& x) const {
        if(!geometryValidator)return true;
        zPointArray points;for(int i=0;i<x.rows();++i)points.emplace_back(float(x(i,0)),float(x(i,1)),float(x(i,2)));
        return geometryValidator(points);
    }
    Vec residualVector(const Mat& x,const Vec& values) const {
        require(validGeometry(x),"Stair clearance or flight layout constraint violated.");
        const int n=int(free.size());
        Vec r=-loadVector();
        for(int j=0;j<int(active.size());++j) {
            int e=active[j], u=edges[2*e], v=edges[2*e+1];
            Eigen::RowVector3d d=x.row(u)-x.row(v);
            require(d.allFinite() && d.norm()>minLength, "Collapsed or nonfinite member.");
            for(int k=0;k<3;++k) {
                if(row[u]>=0) r[row[u]+k*n]+=d[k]*values[j];
                if(row[v]>=0) r[row[v]+k*n]-=d[k]*values[j];
            }
        }
        return r;
    }
    Vec bounded(const Vec& values) const { return values.cwiseMax(lo).cwiseMin(hi); }
    double regularizationWeight() const { return settings.edgeRegularizationStrength*settings.geometryWeight; }
    bool fixedXY() const { return settings.fixXY||settings.planConstraint==zVaultPlanConstraint::FixXY; }
    double planWeight() const { return !fixedXY()&&settings.planConstraint==zVaultPlanConstraint::PenalizeXY?settings.xyPenaltyWeight*settings.geometryWeight:0; }
    double planPenalty(const Mat& x) const {
        double sum=0;for(int i:free){double dx=x(i,0)-reference[i].x,dy=x(i,1)-reference[i].y;sum+=dx*dx+dy*dy;}
        return sum/(free.size()*planReferenceLength*planReferenceLength);
    }
    // RMS length deficit relative to a fixed initial scale. The corresponding
    // negative gradient extends short edges AND shortens the longest edges.
    double edgePenalty(const Mat& x, Mat* correction=nullptr, double* ratio=nullptr) const {
        Vec lengths(active.size());
        for(int j=0;j<int(active.size());++j){int e=active[j];lengths[j]=(x.row(edges[2*e])-x.row(edges[2*e+1])).norm();}
        double longest=lengths.maxCoeff(), shortest=lengths.minCoeff();
        if(ratio)*ratio=longest>minLength?shortest/longest:0;
        if(correction)*correction=Mat::Zero(x.rows(),3);
        double total=0,sum=0;int ties=0;
        for(int j=0;j<lengths.size();++j)if(lengths[j]==longest)++ties;
        for(int j=0;j<lengths.size();++j){double deficit=std::max(0.0,settings.minMaxEdgeRatio*longest-lengths[j]);total+=deficit;sum+=deficit*deficit;}
        if(correction)for(int j=0;j<lengths.size();++j){
            if(lengths[j]<=minLength)continue;
            int e=active[j],u=edges[2*e],v=edges[2*e+1];
            double amount=std::max(0.0,settings.minMaxEdgeRatio*longest-lengths[j]);
            if(lengths[j]==longest)amount-=settings.minMaxEdgeRatio*total/ties;
            Eigen::RowVector3d f=(amount/(lengths[j]*lengths.size()))*(x.row(u)-x.row(v));
            if(row[u]>=0)correction->row(u)+=f;if(row[v]>=0)correction->row(v)-=f;
        }
        return sum/(lengths.size()*edgeReferenceLength*edgeReferenceLength);
    }
    double merit(const Mat& x,double residual) const {
        return residual*residual/std::max(1e-28,loadVector().squaredNorm())+regularizationWeight()*edgePenalty(x)+planWeight()*planPenalty(x);
    }
    Vec leastSquares(const Mat& x,const Vec& rhs) const {
        const int n=int(free.size()), m=int(active.size());
        if(rhs.squaredNorm()<1e-30) return Vec::Zero(m);
        std::vector<Eigen::Triplet<double>> terms;
        terms.reserve(size_t(m)*6);
        for(int j=0;j<m;++j) {
            int e=active[j], u=edges[2*e], v=edges[2*e+1];
            Eigen::RowVector3d d=x.row(u)-x.row(v);
            require(d.allFinite() && d.norm()>minLength, "Collapsed or nonfinite member.");
            for(int k=0;k<3;++k) {
                if(row[u]>=0) terms.emplace_back(row[u]+k*n,j,d[k]);
                if(row[v]>=0) terms.emplace_back(row[v]+k*n,j,-d[k]);
            }
        }
        Eigen::SparseMatrix<double> a(3*n,m);
        a.setFromTriplets(terms.begin(),terms.end());
        if(largeSystem()) {
            // Iterate on sparse A and A^T; never allocate the dense equilibrium
            // matrix or form/factor A^T A for refined meshes.
            Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
            solver.setMaxIterations(1000);solver.setTolerance(1e-8);solver.compute(a);
            Vec solution=solver.solve(rhs);
            require(solution.allFinite(),"Sparse force density fit produced nonfinite values.");
            return solution;
        }
        Eigen::SparseMatrix<double> normal=a.transpose()*a;
        for(int i=0;i<m;++i) normal.coeffRef(i,i)+=0.0;
        normal.makeCompressed();
        Vec b=a.transpose()*rhs, diagonal=normal.diagonal();
        double scale=1.0;
        if(diagonal.size()>0) scale=std::max(scale,diagonal.cwiseAbs().maxCoeff());
        double damping=scale*1e-10;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(normal);
        if(solver.info()==Eigen::Success) for(int attempt=0;attempt<4;++attempt,damping*=10) {
            Eigen::SparseMatrix<double> system=normal;
            for(int i=0;i<m;++i) system.coeffRef(i,i)+=damping*std::max(1.0,diagonal[i]);
            solver.factorize(system);
            if(solver.info()!=Eigen::Success) continue;
            Vec solution=solver.solve(b);
            if(solver.info()==Eigen::Success && solution.allFinite()) return solution;
        }
        Mat dense=matrix(x);
        return dense.colPivHouseholderQr().solve(rhs);
    }
    void gradient(const Mat& x,const Vec& values,Mat& dx,Vec& dq) const {
        Vec r=residualVector(x,values);
        Vec signedSum=Vec::Zero(free.size()), absoluteSum=signedSum;
        for(int j=0;j<int(active.size());++j) {
            int e=active[j], u=edges[2*e], v=edges[2*e+1];
            double f=(x.row(u)-x.row(v)).norm()*values[j];
            for(int id:{u,v}) if(row[id]>=0) { signedSum[row[id]]+=f; absoluteSum[row[id]]+=std::abs(f); }
        }
        dx=Mat::Zero(x.rows(),3);
        for(int i=0;i<int(free.size());++i) {
            // Eq.14 ratio, not sign(sum). Zero-force nodes have zero geometric correction.
            double k=absoluteSum[i]>1e-14?signedSum[i]/absoluteSum[i]:0;
            for(int c=fixedXY()?2:0;c<3;++c) dx(free[i],c)=-k*settings.geometryWeight*r[i+c*free.size()];
        }
        if(regularizationWeight()>0){
            Mat correction;edgePenalty(x,&correction);dx+=regularizationWeight()*correction;
            if(fixedXY())dx.leftCols(2).setZero();
        }
        if(planWeight()>0)for(int i:free){
            dx(i,0)-=planWeight()*(x(i,0)-reference[i].x)/free.size();
            dx(i,1)-=planWeight()*(x(i,1)-reference[i].y)/free.size();
        }
        dq=leastSquares(x,-(1-settings.geometryWeight)*r);
        for(int j=0;j<dq.size();++j) if((values[j]<=lo[j] && dq[j]<0)||(values[j]>=hi[j] && dq[j]>0)) dq[j]=0;
        require(dx.allFinite() && dq.allFinite(), "Nonfinite equilibrium gradient.");
    }
    // Damped Gauss-Newton rescue for a non-descent/cancelled heuristic update.
    // Build the Jacobian of the same normalized merit used for acceptance.
    bool fallback() {
        const int n=int(free.size()), m=int(active.size()), cols=3*n+m;
        const double gx=std::sqrt(settings.geometryWeight), gq=std::sqrt(1-settings.geometryWeight);
        const double loadScale=std::sqrt(std::max(1e-28,loadVector().squaredNorm()));
        std::vector<Eigen::Triplet<double>> terms;
        Vec residual=Vec::Zero(5*n+m);
        residual.head(3*n)=residualVector(X,q)/loadScale;
        auto geometryTerm=[&](int equation,int vertex,int c,double value){
            if(row[vertex]>=0 && (!fixedXY()||c==2) && gx>0)
                terms.emplace_back(equation,c*n+row[vertex],value*gx);
        };
        Vec lengths(m);for(int j=0;j<m;++j){int e=active[j];lengths[j]=(X.row(edges[2*e])-X.row(edges[2*e+1])).norm();}
        double longest=lengths.maxCoeff();int ties=0;for(int j=0;j<m;++j)if(lengths[j]==longest)++ties;
        for(int j=0;j<m;++j){
            int e=active[j],u=edges[2*e],v=edges[2*e+1];
            for(int c=0;c<3;++c)for(int endpoint:{u,v})if(row[endpoint]>=0){
                int equation=c*n+row[endpoint];double sign=endpoint==u?1:-1;
                geometryTerm(equation,u,c,sign*q[j]/loadScale);
                geometryTerm(equation,v,c,-sign*q[j]/loadScale);
                if(gq>0 && lo[j]<hi[j])terms.emplace_back(equation,3*n+j,sign*(X(u,c)-X(v,c))*gq/loadScale);
            }
        }
        double pw=std::sqrt(planWeight()/n)/planReferenceLength;
        for(int i=0;i<n;++i)for(int c=0;c<2;++c){
            int equation=3*n+c*n+i;
            residual[equation]=pw*(X(free[i],c)-(c==0?reference[free[i]].x:reference[free[i]].y));
            geometryTerm(equation,free[i],c,pw);
        }
        double ew=std::sqrt(regularizationWeight()/m)/edgeReferenceLength;
        auto lengthTerm=[&](int equation,int j,double coefficient){
            int e=active[j],u=edges[2*e],v=edges[2*e+1];
            for(int c=0;c<3;++c){double value=coefficient*(X(u,c)-X(v,c))/lengths[j];
                geometryTerm(equation,u,c,value);geometryTerm(equation,v,c,-value);}
        };
        if(ew>0)for(int j=0;j<m;++j){
            double deficit=settings.minMaxEdgeRatio*longest-lengths[j];if(deficit<=0)continue;
            int equation=5*n+j;residual[equation]=ew*deficit;lengthTerm(equation,j,-ew);
            for(int k=0;k<m;++k)if(lengths[k]==longest)lengthTerm(equation,k,ew*settings.minMaxEdgeRatio/ties);
        }
        Eigen::SparseMatrix<double> jac(residual.size(),cols);jac.setFromTriplets(terms.begin(),terms.end());
        Eigen::SparseMatrix<double> normal=jac.transpose()*jac;
        for(int i=0;i<cols;++i) normal.coeffRef(i,i)+=0.0;
        normal.makeCompressed();
        Vec rhs=-jac.transpose()*residual, diagonal=normal.diagonal();
        double damping=1e-4, before=merit(X,report.residualNorm);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        if(!largeSystem()){solver.analyzePattern(normal);if(solver.info()!=Eigen::Success)return false;}
        for(int attempt=0;attempt<12;++attempt,damping*=10){
            Eigen::SparseMatrix<double> system=normal;
            for(int i=0;i<cols;++i)system.coeffRef(i,i)+=damping*std::max(1.0,diagonal[i]);
            Vec delta;
            if(largeSystem()) {
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Lower|Eigen::Upper> iterative;
                iterative.setMaxIterations(1000);iterative.setTolerance(1e-8);iterative.compute(system);delta=iterative.solve(rhs);
                if(!delta.allFinite())continue;
            }else{
                solver.factorize(system);if(solver.info()!=Eigen::Success)continue;
                delta=solver.solve(rhs);if(solver.info()!=Eigen::Success||!delta.allFinite())continue;
            }
            for(double alpha=1;alpha>=1.0/1024;alpha*=0.5){
                Mat trial=X;for(int i=0;i<n;++i)for(int c=fixedXY()?2:0;c<3;++c)trial(free[i],c)+=alpha*gx*delta[c*n+i];
                Vec qt=bounded(q+alpha*gq*delta.tail(m));double norm;
                try{norm=residualVector(trial,qt).norm();}catch(const std::invalid_argument&){continue;}
                if(std::isfinite(norm)&&merit(trial,norm)<before){
                    history.push_back(state());if(history.size()>256)history.pop_front();
                    X=trial;q=qt;report.iteration++;report.status="running";h=settings.initialStep;sync();return true;
                }
            }
        }
        return false;
    }
    State state() const { return {X,q,report,h}; }
    void restore(const State& s) { X=s.X;q=s.q;report=s.report;h=s.h; sync(); }
    void sync() {
        diagrams.report.status = "stale";
        zPointArray positions; positions.reserve(original.size());
        for(int i=0;i<X.rows();++i) positions.emplace_back(float(X(i,0)),float(X(i,1)),float(X(i,2)));
        zFnMesh fn(mesh); fn.setVertexPositions(positions); fn.computeMeshNormals();
        densities.assign(edges.size()/2,0); forces=densities;
        Mat balance=-P;
        for(int j=0;j<int(active.size());++j) {
            int e=active[j],u=edges[2*e],v=edges[2*e+1];
            Eigen::RowVector3d f=q[j]*(X.row(u)-X.row(v));
            balance.row(u)+=f; balance.row(v)-=f;
            densities[e]=q[j]; forces[e]=q[j]*(X.row(u)-X.row(v)).norm();
        }
        residuals.assign(original.size(),zVector()); reactions=residuals;
        double sum=0;report.maxResidual=0;report.maxDisplacement=0;
        for(int i=0;i<X.rows();++i) {
            zVector value(float(balance(i,0)),float(balance(i,1)),float(balance(i,2)));
            if(row[i]<0) reactions[i]=value;
            else { residuals[i]=value*-1; sum+=balance.row(i).squaredNorm(); report.maxResidual=std::max(report.maxResidual,balance.row(i).norm()); }
            Eigen::RowVector3d origin(original[i].x,original[i].y,original[i].z);
            report.maxDisplacement=std::max(report.maxDisplacement,(X.row(i)-origin).norm());
        }
        report.residualNorm=std::sqrt(sum);report.stepSize=h;
        report.rmsXYDisplacement=std::sqrt(planPenalty(X))*planReferenceLength;report.maxXYDisplacement=0;
        for(int i:free)report.maxXYDisplacement=std::max(report.maxXYDisplacement,std::hypot(X(i,0)-reference[i].x,X(i,1)-reference[i].y));
        report.edgeRegularizationError=std::sqrt(edgePenalty(X,nullptr,&report.minMaxEdgeRatio));
        report.converged=report.residualNorm<=std::max(settings.absoluteTolerance,settings.relativeTolerance*report.initialResidual);
        if(regularizationWeight()>0)report.converged=report.converged&&report.edgeRegularizationError<=1e-4;
        // A soft target is not a zero-displacement constraint. Require stationary
        // combined updates as well as physical equilibrium before early stopping.
        if(report.converged&&planWeight()>0){Mat dx;Vec dq;gradient(X,q,dx,dq);
            report.converged=dx.norm()/(planReferenceLength*std::sqrt(double(free.size())))<=settings.integrationTolerance;}
        if(report.converged) report.status="converged";
    }
};
zTsVault::zTsVault():impl(new Impl) {}
zTsVault::~zTsVault()=default;
void zTsVault::clear() { impl.reset(new Impl); }
void zTsVault::setInputMesh(zObjectMesh& mesh) {
    auto next=std::make_unique<Impl>(); zFnMesh fn(mesh); zIntArray counts,connects;
    fn.getVertexPositions(next->original); fn.getEdgeData(next->edges); fn.getPolygonData(connects,counts);
    require(!next->original.empty() && !next->edges.empty(),"Assign a nonempty mesh.");
    require(next->original.size()<size_t(std::numeric_limits<int>::max()/8)&&next->edges.size()<size_t(std::numeric_limits<int>::max()/8),"FDM topology exceeds representable solver indices.");
    for(auto& p:next->original) require(finite(p),"Nonfinite mesh coordinate.");
    next->reference=next->original;
    zFnMesh(next->mesh).create(next->original,counts,connects);
    next->masses.assign(next->original.size(),1);next->loads.assign(next->original.size(),zVector());
    impl=std::move(next);
}
void zTsVault::setGeometryValidator(std::function<bool(const zPointArray&)> validator) {
    impl->geometryValidator=std::move(validator);
}
void zTsVault::setReferencePositions(const zPointArray& positions){
    require(!positions.empty()&&positions.size()==impl->original.size(),"Reference positions must match input vertices.");
    for(const auto& p:positions)require(finite(p),"Nonfinite reference position.");
    impl->reference=positions;impl->invalidate();
}
void zTsVault::setSupports(const zIntArray& ids) { impl->ids(ids,impl->original.size());impl->supports=ids;impl->invalidate(); }
void zTsVault::setTensionEdges(const zIntArray& ids) { impl->ids(ids,impl->edges.size()/2);impl->tension=ids;impl->invalidate(); }
void zTsVault::setVertexMass(double mass) { setVertexMasses(zDoubleArray(impl->original.size(),mass)); }
void zTsVault::setVertexMasses(const zDoubleArray& masses) {
    require(masses.size()==impl->original.size(),"Mass count must match vertices.");
    for(double m:masses) require(std::isfinite(m)&&m>=0,"Vertex mass must be finite and nonnegative.");
    impl->masses=masses;impl->invalidate();
}
void zTsVault::setGravity(const zVector& v) { require(finite(v),"Nonfinite gravity.");impl->gravity=v;impl->invalidate(); }
void zTsVault::setVertexLoads(const zVectorArray& loads) {
    require(loads.size()==impl->original.size(),"Load count must match vertices.");
    for(auto& v:loads) require(finite(v),"Nonfinite load.");impl->loads=loads;impl->invalidate();
}
void zTsVault::setForceDensityBounds(const zDoubleArray& lower,const zDoubleArray& upper) {
    require(lower.size()==upper.size() && (lower.empty() || lower.size()==impl->edges.size()/2),"Bounds must match mesh edges.");
    for(size_t i=0;i<lower.size();++i) require(std::isfinite(lower[i])&&std::isfinite(upper[i])&&lower[i]<=upper[i],"Invalid force density interval.");
    impl->lower=lower;impl->upper=upper;impl->invalidate();
}
void zTsVault::setBestFitSettings(const zVaultBestFitSettings& s) {
    require(std::isfinite(s.qMin)&&std::isfinite(s.qMax)&&s.qMin>=0&&s.qMax>=s.qMin&&s.qMax>0,"Require 0 <= qMin <= qMax and qMax > 0.");
    require(std::isfinite(s.geometryWeight)&&s.geometryWeight>=0&&s.geometryWeight<=1,"Geometry weight must be in [0,1].");
    require(int(s.planConstraint)>=0&&int(s.planConstraint)<=2,"Invalid plan constraint.");
    require(std::isfinite(s.xyPenaltyWeight)&&s.xyPenaltyWeight>=0,"XY penalty weight must be nonnegative.");
    require(std::isfinite(s.minMaxEdgeRatio)&&s.minMaxEdgeRatio>=0&&s.minMaxEdgeRatio<=1,"Min/max edge ratio must be in [0,1].");
    require(std::isfinite(s.edgeRegularizationStrength)&&s.edgeRegularizationStrength>=0,"Edge regularization strength must be nonnegative.");
    require(std::isfinite(s.absoluteTolerance)&&s.absoluteTolerance>0&&std::isfinite(s.relativeTolerance)&&s.relativeTolerance>=0&&s.relativeTolerance<1,"Invalid residual tolerance.");
    require(std::isfinite(s.initialStep)&&s.initialStep>0&&std::isfinite(s.integrationTolerance)&&s.integrationTolerance>0&&s.maxIterations>0,"Invalid integration settings.");
    impl->settings=s;impl->invalidate();
}
zVaultSolveReport zTsVault::computeBestFitForceDensities() {
    auto& s=*impl; s.invalidate(); int n=int(s.original.size());require(n>0,"Assign mesh first.");
    s.row.assign(n,-1);s.free.clear();s.active.clear();
    std::vector<bool> fixed(n,false),tension(s.edges.size()/2,false);
    for(int i:s.supports) fixed[i]=true;for(int e:s.tension) tension[e]=true;
    for(int i=0;i<n;++i) if(!fixed[i]) { s.row[i]=int(s.free.size());s.free.push_back(i); }
    require(!s.free.empty(),"At least one free vertex is required.");
    for(int e=0;e<int(s.edges.size()/2);++e) if(!fixed[s.edges[2*e]]||!fixed[s.edges[2*e+1]]) s.active.push_back(e);
    require(!s.active.empty(),"At least one active edge is required.");
    s.X.resize(n,3);s.P.resize(n,3);
    for(int i=0;i<n;++i) for(int c=0;c<3;++c) {
        s.X(i,c)=c==0?s.original[i].x:c==1?s.original[i].y:s.original[i].z;
        auto load=s.loads[i]+s.gravity*float(s.masses[i]);
        require(finite(load),"Nonfinite mass-generated load.");s.P(i,c)=c==0?load.x:c==1?load.y:load.z;
    }
    s.minLength=std::max(1e-12,(s.X.colwise().maxCoeff()-s.X.colwise().minCoeff()).norm()*1e-10);
    s.edgeReferenceLength=s.minLength;
    for(int e:s.active)s.edgeReferenceLength=std::max(s.edgeReferenceLength,(s.X.row(s.edges[2*e])-s.X.row(s.edges[2*e+1])).norm());
    s.planReferenceLength=1e-12;
    for(int e=0;e<int(s.edges.size()/2);++e)s.planReferenceLength=std::max(s.planReferenceLength,double((s.reference[s.edges[2*e]]-s.reference[s.edges[2*e+1]]).length()));
    s.lo.resize(s.active.size());s.hi.resize(s.active.size());
    for(int j=0;j<int(s.active.size());++j) {
        int e=s.active[j];s.lo[j]=tension[e]?s.settings.qMin:-s.settings.qMax;
        s.hi[j]=tension[e]?s.settings.qMax:s.settings.compressionOnly?-s.settings.qMin:s.settings.qMax;
        if(!s.lower.empty()) {s.lo[j]=std::max(s.lo[j],s.lower[e]);s.hi[j]=std::min(s.hi[j],s.upper[e]);}
        require(s.lo[j]<=s.hi[j],"Force density bounds conflict with edge sign assignment.");
    }
    require(s.loadVector().norm()>1e-14,"Nonzero free-vertex loads required; prestress-only solving is not supported.");
    if(s.largeSystem()) {
        s.q=s.bounded(s.leastSquares(s.X,s.loadVector()));
        s.report.rank=-1; // Iterative least squares does not measure numerical rank.
    }else{
        Mat a=s.matrix(s.X);auto qr=a.completeOrthogonalDecomposition();s.q=s.bounded(qr.solve(s.loadVector()));
        s.report.rank=int(qr.rank());
    }
    require(s.q.allFinite(),"Force density fit failed.");
    s.report.activeEdges=int(s.active.size());s.report.freeVertices=int(s.free.size());
    s.report.status="fitted";s.h=s.settings.initialStep;s.sync();
    s.report.initialResidual=s.report.residualNorm;s.sync();s.initial=s.state();s.initialized=true;
    return s.report;
}
zVaultSolveReport zTsVault::stepBestFit(int iterations) {
    auto& s=*impl;require(s.initialized,"Fit force densities before stepping.");require(iterations>=1&&iterations<=10000,"Invalid iteration count.");
    for(int iter=0;iter<iterations && !s.report.converged;++iter) {
        if(s.report.status=="stalled") break;
        if(s.report.iteration>=s.settings.maxIterations) {s.report.status="iteration-limit";break;}
        auto before=s.state();bool accepted=false;
        // Adaptive explicit RK2 (Heun), with Euler error estimate. Regularized
        // runs use a combined merit; otherwise require physical residual decrease.
        // All stage densities are projected. Numerical departures are documented.
        for(int attempt=0;attempt<24;++attempt) {
            Mat dx1,dx2;Vec dq1,dq2;
            s.gradient(s.X,s.q,dx1,dq1);
            Mat trial=s.X+s.h*dx1;Vec qt=s.bounded(s.q+s.h*dq1);
            try { s.gradient(trial,qt,dx2,dq2); } catch(const std::invalid_argument&) {s.h*=0.5;continue;}
            Mat x=s.X+s.h*0.5*(dx1+dx2);Vec q=s.bounded(s.q+s.h*0.5*(dq1+dq2));
            double error=std::max((x-trial).norm()/(1+x.norm()),(q-qt).norm()/(1+q.norm()));
            double norm;
            try {norm=s.residualVector(x,q).norm();} catch(const std::invalid_argument&) {s.h*=0.5;continue;}
            bool improves=s.regularizationWeight()>0||s.planWeight()>0?s.merit(x,norm)<s.merit(s.X,s.report.residualNorm):norm<s.report.residualNorm;
            if(std::isfinite(norm)&&error<=s.settings.integrationTolerance&&improves) {
                s.history.push_back(before);if(s.history.size()>256) s.history.pop_front();
                s.X=x;s.q=q;s.report.iteration++;s.report.status="running";
                s.h=std::min(s.h*(error<s.settings.integrationTolerance*0.25?1.5:1.0),s.settings.initialStep*16);
                s.sync();accepted=true;break;
            }
            s.h*=0.5;if(s.h<1e-12) break;
        }
        if(!accepted && !s.fallback()) {s.report.status="stalled";break;}
    }
    return s.report;
}
zVaultSolveReport zTsVault::solveBestFit(int count) {return stepBestFit(count);}
zVaultSolveReport zTsVault::snapFixedDensities(double geometryWeight) {
    auto& s=*impl;require(s.initialized,"Fit force densities before snapping.");
    require(std::isfinite(geometryWeight)&&geometryWeight>=0&&geometryWeight<=1,"Snap geometry weight must be in [0,1].");
    const int n=int(s.free.size());std::vector<Eigen::Triplet<double>> terms;terms.reserve(s.active.size()*4);
    Mat rhs(n,3);for(int i=0;i<n;++i)rhs.row(i)=s.P.row(s.free[i]);
    for(int j=0;j<int(s.active.size());++j){int e=s.active[j],u=s.edges[2*e],v=s.edges[2*e+1];double q=s.q[j];
        for(auto endpoints:{std::pair<int,int>{u,v},{v,u}}){int a=s.row[endpoints.first],b=s.row[endpoints.second];if(a<0)continue;
            terms.emplace_back(a,a,q);if(b>=0)terms.emplace_back(a,b,-q);else rhs.row(a)+=q*s.X.row(endpoints.second);
        }
    }
    Eigen::SparseMatrix<double> matrix(n,n);matrix.setFromTriplets(terms.begin(),terms.end());
    Mat reference(n,3);for(int i=0;i<n;++i){auto p=s.reference[s.free[i]];reference.row(i)<<p.x,p.y,p.z;}
    Mat z=reference;
    if(geometryWeight<1) {
    Eigen::SparseMatrix<double> system=matrix;Mat target=rhs;
    if(geometryWeight>0){
        system=matrix.transpose()*matrix;
        double scale=std::max(1e-12,system.diagonal().mean());double lambda=scale*geometryWeight/(1-geometryWeight);
        for(int i=0;i<n;++i)system.coeffRef(i,i)+=lambda;
        target=matrix.transpose()*rhs+lambda*reference;
    }
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;solver.compute(system);
    require(solver.info()==Eigen::Success,"Fixed-density system is singular; check positive density magnitude and supports.");
    double scale=std::max(1.0,solver.vectorD().cwiseAbs().maxCoeff());
    require(solver.vectorD().cwiseAbs().minCoeff()>scale*1e-12,"Fixed-density system is singular; check positive density magnitude and supports.");
    z=solver.solve(target);require(solver.info()==Eigen::Success&&z.allFinite()&&(system*z-target).norm()<=1e-7*(1+target.norm()),"Direct height solve failed numerical verification.");
    }
    Mat candidate=s.X;for(int i=0;i<n;++i)if(s.fixedXY())candidate(s.free[i],2)=z(i,2);else candidate.row(s.free[i])=z.row(i);
    bool shortened=false;
    if(!s.validGeometry(candidate)) {
        const Mat proposed=candidate;bool feasible=false;
        for(int attempt=1;attempt<=30;++attempt){candidate=s.X+std::ldexp(1.0,-attempt)*(proposed-s.X);if(s.validGeometry(candidate)){feasible=true;shortened=true;break;}}
        if(!feasible){s.report.status="clearance-blocked";return s.report;}
    }
    s.residualVector(candidate,s.q); // Validate members before publishing.
    s.history.push_back(s.state());if(s.history.size()>256)s.history.pop_front();s.X=std::move(candidate);
    s.report.iteration++;s.report.status=shortened?"direct-constrained":s.fixedXY()?"direct-height":"direct-equilibrium";s.sync();return s.report;
}
void zTsVault::resetBestFit() {require(impl->initialized,"Fit before reset.");impl->restore(impl->initial);impl->history.clear();}
bool zTsVault::previousBestFit() {require(impl->initialized,"Fit before rewind.");if(impl->history.empty())return false;auto state=impl->history.back();impl->history.pop_back();impl->restore(state);return true;}
zObjectMesh& zTsVault::resultMesh(){return impl->mesh;}
const zDoubleArray& zTsVault::forceDensities()const{return impl->densities;}
const zDoubleArray& zTsVault::axialForces()const{return impl->forces;}
const zVectorArray& zTsVault::residualForces()const{return impl->residuals;}
const zVectorArray& zTsVault::supportReactions()const{return impl->reactions;}
const zIntArray& zTsVault::activeEdgeIds()const{return impl->active;}
const zIntArray& zTsVault::edgeConnects()const{return impl->edges;}
const zVaultSolveReport& zTsVault::report()const{return impl->report;}
zVaultDiagramReport zTsVault::createFormFromResult() {
    require(impl->initialized,"Fit force densities before creating diagrams.");
    zPointArray p; zIntArray counts, connects; zFnMesh fn(impl->mesh);
    fn.getVertexPositions(p); fn.getPolygonData(connects,counts);
    for(auto& v:p)v.z=0;
    zFnMesh(impl->diagrams.form).create(p,counts,connects);
    impl->diagrams.report = {}; impl->diagrams.report.status="form-created";
    return impl->diagrams.report;
}
zVaultDiagramReport zTsVault::createForceFromForm(const zVaultDiagramSettings& settings) {
    require(impl->initialized,"Fit force densities before creating diagrams.");
    for(int i:impl->free) require(impl->P.row(i).head(2).norm()<=1e-10,"Stage 1 diagrams require vertical loads at free vertices.");
    impl->diagrams=buildVaultDiagrams(impl->mesh,impl->supports,impl->edges,impl->densities,settings);
    return impl->diagrams.report;
}
zVaultDiagramReport zTsVault::updateForceFromDensities(){return createForceFromForm(impl->diagrams.settings);}
const zObjectMesh& zTsVault::formMesh()const{return impl->diagrams.form;}
const zObjectMesh& zTsVault::forceMesh()const{return impl->diagrams.force;}
const zIntArray& zTsVault::formHalfEdgeToForceHalfEdge()const{return impl->diagrams.formToForce;}
const zIntArray& zTsVault::forceHalfEdgeToFormHalfEdge()const{return impl->diagrams.forceToForm;}
const zVaultDiagramReport& zTsVault::diagramReport()const{return impl->diagrams.report;}
}
