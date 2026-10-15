#include <zspace/zToolsets/fabrication/zTsCTubeBranching.h>
#include <Eigen/Dense>
#include "../../../third_party/LBFGSpp/include/LBFGSB.h"
#include <map>
#include <tuple>
#include <cmath>
#include <stdexcept>
namespace zSpace {
namespace {
using V=Eigen::Vector3d;
V vec(const zPoint& p){return V(p.x,p.y,p.z);}
zPoint point(const V& p){return zPoint(float(p.x()),float(p.y()),float(p.z()));}
V unit(const V& p){if(!p.allFinite()||p.norm()<1e-10)throw std::invalid_argument("Degenerate central optimisation frame.");return p.normalized();}
struct Panel {int a,b;V origin,step;};
}
void zTsCTubeBranching::optimiseCentral(const std::vector<zCTubeBranchInput>& input,const zCTubeBranchSettings& settings) {
    // Coplanarity is symmetric between the two rails. For constrained boundary
    // optimisation, use the central rail as fixed and the boundary as unknown.
    auto solvingInput=input;
    if(!settings.fixedDirectrices)for(auto& branch:solvingInput){if(branch.frames.size()!=branch.centralTargets.size())throw std::invalid_argument("Boundary optimisation needs matching frame and seam stations.");for(size_t j=0;j<branch.frames.size();++j)std::swap(branch.frames[j].origin,branch.centralTargets[j]);}
    std::map<std::tuple<float,float,float>,int> ids;
    std::vector<V> initial;std::vector<std::vector<int>> branchIds;
    std::vector<Panel> panels;double scale=0;int count=0;
    for(const auto& branch:solvingInput){
        if(branch.frames.size()<2||branch.frames.size()>512||branch.frames.size()!=branch.centralTargets.size())throw std::invalid_argument("Central optimisation needs matching frame and seam stations.");
        branchIds.emplace_back();
        for(size_t j=0;j<branch.frames.size();++j){const auto p=branch.centralTargets[j];const auto key=std::make_tuple(p.x,p.y,p.z);auto found=ids.find(key);if(found==ids.end()){const int id=int(initial.size());ids.emplace(key,id);initial.push_back(vec(p));found=ids.find(key);}branchIds.back().push_back(found->second);scale+=(vec(p)-vec(branch.frames[j].origin)).norm();++count;
            if(j){const V origin=vec(branch.frames[j-1].origin);panels.push_back({branchIds.back()[j-1],branchIds.back()[j],origin,vec(branch.frames[j].origin)-origin});}
        }
    }
    scale/=count;if(scale<1e-8)throw std::invalid_argument("Collapsed central flange.");
    const int P=3*int(initial.size());Eigen::VectorXd q(P);for(size_t i=0;i<initial.size();++i)q.segment<3>(3*i)=initial[i];
    // Shared seam vertices are single unknowns. Scalar triple products enforce
    // panel coplanarity with the fixed boundary edges.
    std::vector<double> multipliers(panels.size(),0);double penalty=10.0/panels.size();
    const double movementWeight=1e-3/(initial.size()*scale*scale),smoothWeight=1e-4/(panels.size()*scale*scale);
    auto evaluate=[&](const Eigen::VectorXd& controls,Eigen::VectorXd& grad,bool gradient){
        if(gradient)grad=Eigen::VectorXd::Zero(P);double energy=0;
        for(size_t i=0;i<initial.size();++i){const V d=controls.segment<3>(3*i)-initial[i];energy+=movementWeight*d.squaredNorm();if(gradient)grad.segment<3>(3*i)+=2*movementWeight*d;}
        for(size_t i=0;i<panels.size();++i){const auto& p=panels[i];const V a=controls.segment<3>(3*p.a)-p.origin,b=controls.segment<3>(3*p.b)-p.origin;
            const double denominator=std::max(1e-10,p.step.norm()*scale),r=p.step.cross(a).dot(b)/denominator;
            energy+=multipliers[i]*r+.5*penalty*r*r;
            if(gradient){const double factor=(multipliers[i]+penalty*r)/denominator;grad.segment<3>(3*p.a)+=factor*b.cross(p.step);grad.segment<3>(3*p.b)+=factor*p.step.cross(a);}
            const V smooth=(controls.segment<3>(3*p.b)-initial[p.b])-(controls.segment<3>(3*p.a)-initial[p.a]);energy+=smoothWeight*smooth.squaredNorm();if(gradient){grad.segment<3>(3*p.a)-=2*smoothWeight*smooth;grad.segment<3>(3*p.b)+=2*smoothWeight*smooth;}
        }
        if(!std::isfinite(energy))throw std::invalid_argument("Central optimisation diverged.");return energy;
    };
    Eigen::VectorXd grad;zCTubeBranchReport report;const double first=evaluate(q,grad,true);report.initialResidual=std::sqrt(std::max(0.0,first));
    for(int k=0;k<std::min(P,12);++k){const int index=k*(P-1)/std::max(1,std::min(P,12)-1);Eigen::VectorXd plus=q,minus=q,dummy;plus[index]+=1e-6;minus[index]-=1e-6;const double fd=(evaluate(plus,dummy,false)-evaluate(minus,dummy,false))/2e-6;report.gradientCheckError=std::max(report.gradientCheckError,std::abs(fd-grad[index])/std::max(1.0,std::abs(fd)));}
    Eigen::VectorXd lower=q.array()-10*scale,upper=q.array()+10*scale;
    if(settings.fixXY)for(int i=0;i<P;i+=3){lower[i]=upper[i]=q[i];lower[i+1]=upper[i+1]=q[i+1];}
    // Analytic Gauss-Newton projection gives the constrained L-BFGS solves a
    // feasible starting graph, rather than driving an ill-conditioned penalty
    // objective from the imported nonplanar ribbon.
    auto projectPlanarity=[&](Eigen::VectorXd& controls){
        for(int iteration=0;iteration<80;++iteration){
            Eigen::VectorXd residual(panels.size());Eigen::MatrixXd J=Eigen::MatrixXd::Zero(panels.size(),P);double maximum=0;
            for(size_t i=0;i<panels.size();++i){const auto& p=panels[i];const V a=controls.segment<3>(3*p.a)-p.origin,b=controls.segment<3>(3*p.b)-p.origin,n=p.step.cross(a);const double denominator=std::max(1e-10,p.step.norm()*scale);residual[i]=n.dot(b)/denominator;J.block<1,3>(i,3*p.a)=b.cross(p.step).transpose()/denominator;J.block<1,3>(i,3*p.b)=n.transpose()/denominator;maximum=std::max(maximum,std::abs(n.dot(b))/std::max(n.norm(),1e-10));}
            if(settings.fixXY)for(int i=0;i<P;i+=3){J.col(i).setZero();J.col(i+1).setZero();}
            if(maximum<settings.tolerance*.01)break;
            Eigen::MatrixXd system=J*J.transpose();system.diagonal().array()+=1e-10;const Eigen::VectorXd step=J.transpose()*system.ldlt().solve(-residual);
            bool accepted=false;
            for(double fraction=1;fraction>1e-6;fraction*=.5){const Eigen::VectorXd candidate=(controls+fraction*step).cwiseMax(lower).cwiseMin(upper);double error=0;for(const auto& p:panels){const V a=candidate.segment<3>(3*p.a)-p.origin,b=candidate.segment<3>(3*p.b)-p.origin;const double r=p.step.cross(a).dot(b)/std::max(1e-10,p.step.norm()*scale);error+=r*r;}if(error<residual.squaredNorm()){controls=candidate;accepted=true;break;}}
            if(!accepted)break;
        }
    };
    projectPlanarity(q);
    for(int outer=0;outer<10;++outer){
        Eigen::VectorXd best=q;double bestEnergy=evaluate(q,grad,true);
        auto objective=[&](const Eigen::VectorXd& controls,Eigen::VectorXd& g){const double e=evaluate(controls,g,true);if(e<bestEnergy){bestEnergy=e;best=controls;}return e;};
        LBFGSpp::LBFGSBParam<double> options;options.m=16;options.max_iterations=settings.iterations;options.epsilon=1e-10;options.epsilon_rel=1e-10;options.max_linesearch=60;
        LBFGSpp::LBFGSBSolver<double> solver(options);double energy=bestEnergy;try{report.iterations+=solver.minimize(objective,q,energy,lower,upper);}catch(const std::runtime_error&){}q=best;
        report.maxResidual=0;
        for(size_t i=0;i<panels.size();++i){const auto& p=panels[i];const V a=q.segment<3>(3*p.a)-p.origin,b=q.segment<3>(3*p.b)-p.origin,n=p.step.cross(a);const double r=n.dot(b)/std::max(1e-10,p.step.norm()*scale);multipliers[i]+=penalty*r;report.maxResidual=std::max(report.maxResidual,std::abs(n.dot(b))/std::max(n.norm(),1e-10));}
        if(report.maxResidual<=settings.tolerance)break;penalty*=10;
    }
    projectPlanarity(q);
    report.maxResidual=0;for(const auto& p:panels){const V a=q.segment<3>(3*p.a)-p.origin,b=q.segment<3>(3*p.b)-p.origin,n=p.step.cross(a);report.maxResidual=std::max(report.maxResidual,std::abs(n.dot(b))/std::max(n.norm(),1e-10));}
    std::vector<zCTubeBranchResult> output;double movement=0;
    for(size_t i=0;i<initial.size();++i){const double d=(q.segment<3>(3*i)-initial[i]).norm();movement+=d*d;report.maxDisplacement=std::max(report.maxDisplacement,d);}
    for(size_t b=0;b<input.size();++b){const auto& source=input[b];zCTubeBranchResult result;std::vector<V> vertices;
        for(size_t j=0;j<source.frames.size();++j){const V moving=q.segment<3>(3*branchIds[b][j]),O=settings.fixedDirectrices?vec(source.frames[j].origin):moving,C=settings.fixedDirectrices?moving:vec(source.centralTargets[j]),x=unit(C-O),surface=vec(source.frames[j].y),y=unit(surface-x*x.dot(surface)),normal=unit(x.cross(y));
            if(j==0){vertices={O,O-x*settings.outwardLength,C,O-y*settings.stemLength};}
            else {
                const V previousOrigin=vertices[0],step=O-previousOrigin,previousStem=vertices[3]-previousOrigin;
                std::vector<V> next(4);next[0]=O;
                // Coplanar boundary/seam panels also contain the extended
                // flange. Keep its outward length and side at every frame,
                // rather than crossing an apex and reversing the T crossbar.
                next[1]=O-x*settings.outwardLength;
                next[2]=C;
                // The next stem ruling lies in both the current section plane
                // and the plane through the previous stem and boundary edge.
                // Its length is free within that plane, so keep the input stem
                // length rather than shrinking it through a conical apex.
                const V panelNormal=step.cross(previousStem);V direction=panelNormal.cross(normal);
                if(direction.norm()<1e-10)direction=previousStem-normal*normal.dot(previousStem);
                direction=unit(direction);if(direction.dot(-surface)<0)direction=-direction;
                next[3]=O+direction*settings.stemLength;vertices=std::move(next);
            }
            result.frames.push_back({point(O),point(x),point(y),point(normal)});for(const auto& p:vertices)result.sections.push_back(point(p));
        }
        output.push_back(std::move(result));
    }
    report.rmsDisplacement=std::sqrt(movement/initial.size());report.converged=report.maxResidual<=settings.tolerance;m_branches=std::move(output);m_report=report;
}
}
