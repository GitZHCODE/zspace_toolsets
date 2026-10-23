#pragma once
#include <zspace/zToolsets/fabrication/zTsCTubes.h>
#include <zspace/zToolsets/statics/zTsForcePolyhedra.h>
#include <memory>

namespace zSpace {
enum class zForceCTubeMaterial { SolidPolygon, ThinWall };
struct zForceCTubeSettings {
    zForceCTubeMaterial material = zForceCTubeMaterial::ThinWall;
    double thickness = 0.1; // Fixed wall thickness; strip model, junction overlaps neglected.
    double forceScale = 1, referenceStress = 1; // Target material area = force-face area / both factors.
    double areaTolerance = 0.01; // Relative section-area deviation.
    double correspondingAreaTolerance = 1e-4; // Absolute force-face difference.
    double smoothness = 0.001;
    double minScale = 0.05, maxScale = 20;
    int maxIterations = 50;
};
struct zForceCTubeSection {
    double area=0, targetArea=0, relativeAreaError=0;
    zPoint centroid, thrustPoint;
    double eccentricity=0, ex=0, ey=0, kernUtilisation=0;
    double minStress=0, maxStress=0, stressUtilisation=0;
    bool intersectionValid=false, compressionOnly=false, withinStress=false;
    std::string status;
};
struct zForceCTubeReport {
    bool areaResolved=false, forceAreasResolved=false, sourceEquilibrium=false;
    bool sectionChecksPassed=false;
    int iterations=0;
    double initialAreaError=0, maxAreaError=0, maxEccentricity=0;
    double maxQuadWarp=0; // Dimensionless distance of fourth corner from first-triangle plane.
    std::string status="uninitialized", termination;
};
// Frozen-force sizing for one open, sequentially indexed chain. No self-weight feedback.
// Forces are compression-negative. Frames and C-Tube projection construction are retained;
// optimisation varies seed scale and one apex scale per segment, not independent section vertices.
class ZSPACE_TOOLSETS zTsForceCTubes {
public:
    zTsForceCTubes();
    ~zTsForceCTubes();
    zTsForceCTubes(const zTsForceCTubes&)=delete;
    zTsForceCTubes& operator=(const zTsForceCTubes&)=delete;
    void setInput(zObjectGraph& directrix,zObjectGraph& generatrix);
    void setForceDiagram(const zTsForcePolyhedra& force, double forceScale=1);
    void setForceData(zObjectGraph& form,const zDoubleArray& memberFaceAreas,
        const zDoubleArray& axialForces,const zDoubleArray& pairedAreaDifferences,
        bool sourceEquilibrium=false);
    void setSettings(const zForceCTubeSettings& settings);
    zForceCTubeReport compute(bool optimise=false);
    zObjectMesh& tubeMesh();
    const std::vector<zCTubeFrame>& frames() const;
    const std::vector<zForceCTubeSection>& sections() const;
    const zDoubleArray& apexScales() const;
    const zForceCTubeReport& report() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
}
