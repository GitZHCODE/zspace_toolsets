#pragma once
#include <zspace/zToolsets/fabrication/zTsCTubes.h>
namespace zSpace {
struct zCTubeBranchInput {std::vector<zCTubeFrame> frames;zPointArray centralTargets;};
struct zCTubeBranchSettings {double outwardLength=.2,stemLength=.4,tolerance=1e-6;int iterations=400;bool fixedDirectrices=false,fixXY=false;};
struct zCTubeBranchResult {std::vector<zCTubeFrame> frames;zPointArray sections;};
struct zCTubeBranchReport {bool converged=false;int iterations=0;double initialResidual=0,maxResidual=0,maxDisplacement=0,rmsDisplacement=0,gradientCheckError=0,seamError=0;};
class ZSPACE_TOOLSETS zTsCTubeBranching {
public:
    void optimise(const std::vector<zCTubeBranchInput>& input,const zCTubeBranchSettings& settings);
    const std::vector<zCTubeBranchResult>& branches() const {return m_branches;}
    const zCTubeBranchReport& report() const {return m_report;}
private:
    void optimiseCentral(const std::vector<zCTubeBranchInput>& input,const zCTubeBranchSettings& settings);
    std::vector<zCTubeBranchResult> m_branches;zCTubeBranchReport m_report;
};
}
