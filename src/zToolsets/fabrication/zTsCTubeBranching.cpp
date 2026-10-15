#include <zspace/zToolsets/fabrication/zTsCTubeBranching.h>
#include <Eigen/Dense>
#include "../../../third_party/LBFGSpp/include/LBFGSB.h"
#include <cmath>
#include <stdexcept>
namespace zSpace {
namespace {
using V=Eigen::Vector3d;using D=Eigen::Matrix<double,3,Eigen::Dynamic>;
V vec(const zPoint& p){return V(p.x,p.y,p.z);}zPoint point(const V& p){return zPoint(float(p.x()),float(p.y()),float(p.z()));}
V unit(V p){if(!p.allFinite()||p.norm()<1e-10)throw std::invalid_argument("Degenerate branching frame.");return p.normalized();}
}
void zTsCTubeBranching::optimise(const std::vector<zCTubeBranchInput>& input,const zCTubeBranchSettings& settings) {
    if(input.size()!=3||settings.outwardLength<=0||settings.stemLength<=0||settings.iterations<1||settings.tolerance<=0)throw std::invalid_argument("Invalid branching C-Tube settings.");
    if(settings.fixedDirectrices||settings.fixXY){optimiseCentral(input,settings);return;}
    std::vector<zCTubeBranchResult> output;zCTubeBranchReport report;double sumDisplacement=0;int count=0;
    for(const auto& branch:input) {
        const int M=int(branch.frames.size()),P=3*(M-1);
        if(M<2||M>512||branch.centralTargets.size()!=size_t(M))throw std::invalid_argument("Branch optimiser needs 2 to 512 shared seam stations.");
        std::vector<V> C,O,X,Y,N;double scale=0;
        for(int j=0;j<M;++j){C.push_back(vec(branch.centralTargets[j]));O.push_back(vec(branch.frames[j].origin));X.push_back(unit(C.back()-O.back()));Y.push_back(unit(vec(branch.frames[j].y)));N.push_back(unit(X.back().cross(Y.back())));scale+=(C.back()-O.back()).norm();}
        scale/=M;if(scale<1e-8)throw std::invalid_argument("Collapsed branch flange.");
        double extent=0;for(const auto& p:O)extent=std::max(extent,(p-O.front()).norm());
        const double minimumEdge=4e-5*std::max(extent,1e-6);
        std::vector<V> vertices,normals;std::vector<D> derivatives,normalDerivatives;
        auto evaluate=[&](const Eigen::VectorXd& q,Eigen::VectorXd& grad,bool gradient)->double {
            vertices.resize(4*M);normals.resize(M);if(gradient){derivatives.assign(4*M,D::Zero(3,P));normalDerivatives.assign(M,D::Zero(3,P));grad=Eigen::VectorXd::Zero(P);}
            normals[0]=N[0];vertices[0]=O[0];vertices[1]=O[0]-X[0]*settings.outwardLength;vertices[2]=C[0];vertices[3]=O[0]-Y[0]*settings.stemLength;
            for(int j=1;j<M;++j){const int c=3*(j-1);V raw=N[j]+q[c+1]*X[j]+q[c+2]*Y[j];normals[j]=unit(raw);
                if(gradient){Eigen::Matrix3d projection=(Eigen::Matrix3d::Identity()-normals[j]*normals[j].transpose())/raw.norm();normalDerivatives[j].col(c+1)=projection*X[j];normalDerivatives[j].col(c+2)=projection*Y[j];}
                const double alpha=std::exp(q[c]);
                for(int k=0;k<4;++k){V prev=vertices[4*(j-1)+k],r=alpha*(prev-C[j-1]),direction=C[j]+r-prev;double denominator=normals[j].dot(direction),numerator=normals[j].dot(r);
                    if(std::abs(denominator)<1e-8*std::max(1.0,direction.norm()))throw std::invalid_argument("Branch projection parallel to section plane.");
                    double fraction=numerator/denominator;vertices[4*j+k]=C[j]+r-fraction*direction;
                    if(gradient && k==0){D dr=alpha*derivatives[4*(j-1)+k];dr.col(c)+=r;D dd=dr-derivatives[4*(j-1)+k];
                        Eigen::RowVectorXd ds=((r.transpose()*normalDerivatives[j]+normals[j].transpose()*dr)*denominator-(direction.transpose()*normalDerivatives[j]+normals[j].transpose()*dd)*numerator)/(denominator*denominator);
                        derivatives[4*j+k]=dr-fraction*dd-direction*ds;}
                }
                // The seam is a hard construction invariant, not a fitting objective.
                vertices[4*j+2]=C[j];if(gradient)derivatives[4*j+2].setZero();
            }
            // Keep the construction inside the kernel's nondegenerate domain.
            // Extra iterations must not shrink a fabrication edge below its
            // strict unfolding precision while improving the boundary fit.
            for(int j=0;j<M;++j) {
                for(int k=1;k<4;++k)if((vertices[4*j+k]-vertices[4*j]).norm()<minimumEdge)throw std::invalid_argument("Collapsed branch section edge.");
                if(j)for(int k=0;k<4;++k)if((vertices[4*j+k]-vertices[4*(j-1)+k]).norm()<minimumEdge)throw std::invalid_argument("Collapsed branch longitudinal edge.");
            }
            double energy=0;const double fitWeight=1/(M*scale*scale),planeWeight=1e-4/M,smoothWeight=1e-5/M;
            for(int j=0;j<M;++j){V d=vertices[4*j]-O[j],dn=normals[j]-N[j];energy+=fitWeight*d.squaredNorm()+planeWeight*dn.squaredNorm();if(gradient)grad+=2*fitWeight*derivatives[4*j].transpose()*d+2*planeWeight*normalDerivatives[j].transpose()*dn;}
            for(int j=1;j<M-1;++j){V dn=normals[j-1]-2*normals[j]+normals[j+1];energy+=smoothWeight*dn.squaredNorm();if(gradient)grad+=2*smoothWeight*(normalDerivatives[j-1]-2*normalDerivatives[j]+normalDerivatives[j+1]).transpose()*dn;}
            for(int j=1;j<M-2;++j){double a=std::exp(q[3*(j-1)]),b=std::exp(q[3*j]),c=std::exp(q[3*(j+1)]),d=a-2*b+c;energy+=smoothWeight*d*d;if(gradient){grad[3*(j-1)]+=2*smoothWeight*d*a;grad[3*j]-=4*smoothWeight*d*b;grad[3*(j+1)]+=2*smoothWeight*d*c;}}
            if(!std::isfinite(energy)||energy>1e20)throw std::invalid_argument("Branch construction diverged.");return energy;
        };
        Eigen::VectorXd q=Eigen::VectorXd::Zero(P),gradient;double initial=evaluate(q,gradient,true);report.initialResidual=std::max(report.initialResidual,std::sqrt(initial));
        // Check selected components against central differences on the actual construction.
        for(int k=0;k<std::min(P,12);++k){int index=k*(P-1)/std::max(1,std::min(P,12)-1);Eigen::VectorXd plus=q,minus=q,dummy;plus[index]+=1e-6;minus[index]-=1e-6;double fd=(evaluate(plus,dummy,false)-evaluate(minus,dummy,false))/2e-6;report.gradientCheckError=std::max(report.gradientCheckError,std::abs(fd-gradient[index])/std::max(1.0,std::abs(fd)));}
        Eigen::VectorXd best=q;double bestEnergy=initial;
        auto objective=[&](const Eigen::VectorXd& controls,Eigen::VectorXd& grad){try{double e=evaluate(controls,grad,true);if(e<bestEnergy){bestEnergy=e;best=controls;}return e;}catch(const std::invalid_argument&){grad=Eigen::VectorXd::Zero(P);return 1e30;}};
        LBFGSpp::LBFGSBParam<double> options;options.m=12;options.max_iterations=settings.iterations;options.epsilon=1e-8;options.epsilon_rel=1e-8;options.max_linesearch=60;
        LBFGSpp::LBFGSBSolver<double> solver(options);Eigen::VectorXd lower=Eigen::VectorXd::Constant(P,-2),upper=Eigen::VectorXd::Constant(P,2);double energy=initial;
        try{report.iterations+=solver.minimize(objective,q,energy,lower,upper);}catch(const std::runtime_error&){} // Keep the best valid construction after a stalled line search.
        evaluate(best,gradient,true);zCTubeBranchResult result;
        for(int j=0;j<M;++j){V origin=vertices[4*j],x=unit(vertices[4*j+2]-origin),normal=normals[j],y=unit(normal.cross(x));result.frames.push_back({point(origin),point(x),point(y),point(normal)});double distance=(origin-O[j]).norm();report.maxDisplacement=std::max(report.maxDisplacement,distance);sumDisplacement+=distance*distance;++count;for(int k=0;k<4;++k)result.sections.push_back(point(vertices[4*j+k]));}
        output.push_back(std::move(result));report.maxResidual=std::max(report.maxResidual,std::sqrt(bestEnergy));
    }
    report.rmsDisplacement=std::sqrt(sumDisplacement/count);report.converged=report.rmsDisplacement<=settings.tolerance;m_report=report;m_branches=std::move(output);
}
}
