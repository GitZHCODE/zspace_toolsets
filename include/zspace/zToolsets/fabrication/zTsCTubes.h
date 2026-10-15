#pragma once
#include <zspace/interface.h>
#include <zspace/zToolsets/zToolsetsExport.h>

namespace zSpace
{
    enum class zCTubeMode { Generic, ParallelTransport };

    struct zCTubeSettings
    {
        int samples = 64; // Includes the repeated last station for closed curves.
        bool spline = true; // False preserves the input polyline stations.
        double twistRadians = 0.0;
        double scale = 1.0;
        double apexScale = 1.0;
        int fitIterations = 200;
        bool checkFitGradient = false; // Expensive central-difference diagnostic, not used by the solver.
        double fitTolerance = 0.001;
        bool fitPlaneNormals = true;
        double targetSectionWeight = 1.0, targetPlaneWeight = 1.0;
        double apexSmoothnessWeight = 0.0, planeSmoothnessWeight = 0.0;
        zDoubleArray sectionWeights, planeWeights; // Optional per-input multipliers; empty means 1.
        bool preserveOrientation = false; // True keeps the supplied profile plane and in-plane rotation.
    };

    struct zCTubeFitReport {
        bool active=false, converged=false;
        int iterations=0, evaluations=0;
        double gradientNorm=0, gradientCheckError=0;
        std::string termination;
        double initialRms=0, rms=0, maxError=0;
        double initialObjective=0, objective=0, sectionEnergy=0, planeEnergy=0, apexSmoothnessEnergy=0, planeSmoothnessEnergy=0;
        zDoubleArray sectionRms, planeAngleDegrees;
        zIntArray stations;
        zDoubleArray apexScales;
    };
    struct zCTubeFrame { zPoint origin; zVector x, y, normal; };
    struct zCTubeStripArea { double sourceArea, unrolledArea, difference, deviationPercent; };

    // Directrices are chains or loops; sections may be connected planar branched graphs.
    // Each section edge produces a strip. Its centroid or selected anchor is placed on the directrix.
    class ZSPACE_TOOLSETS zTsCTubes
    {
    public:
        void setDirectrix(zObjectGraph& graph);
        void setGeneratrix(zObjectGraph& graph);
        void setGeneratrices(zObjectGraphArray& graphs);
        void setGeneratrixOrigin(int profileIndex, const zPoint& vertex);
        const zCTubeFitReport& fitReport() const { return m_fit; }
        void setSettings(const zCTubeSettings& settings);
        void computeTube(zCTubeMode mode = zCTubeMode::Generic);
        // Fixed mesh-derived section planes; origins must follow the open directrix.
        // The seed profile remains in world coordinates (no refitting or target optimisation).
        void computeTubeWithFrames(const std::vector<zCTubeFrame>& frames,
            const zPointArray& flangeTargets = {}, int flangeVertex = -1);
        void computeTubeWithSections(const std::vector<zCTubeFrame>& frames, const zPointArray& sections);
        void computeStrips();
        void unrollStrips(double gap = 0.1);
        void clearResults();
        const std::vector<zCTubeFrame>& frames() const { return m_frames; }
        const std::vector<zCTubeStripArea>& stripAreas() const { return m_stripAreas; }
        zObjectMesh& tubeMesh() { return m_tube; }
        zObjectGraph& sampledDirectrix() { return m_sampledDirectrix; }
        zObjectMeshArray& stripMeshes() { return m_strips; }
        zObjectMeshArray& unrolledMeshes() { return m_unrolled; }
    private:
        void computeFittedTube();
        std::vector<zPointArray> m_targets;
        zPointArray m_worldProfile;
        zPointArray m_origins;
        std::vector<bool> m_hasOrigin;
        zCTubeFitReport m_fit;
        zPointArray m_curve, m_profile, m_vertices;
        zIntArray m_profileEdges;
        std::vector<std::vector<zIntArray>> m_targetMappings;
        bool m_profileBranched = false;
        zVector m_profileX, m_profileY, m_profileNormal;
        std::vector<zCTubeFrame> m_frames;
        std::vector<zCTubeStripArea> m_stripAreas;
        bool m_closed = false, m_profileClosed = true;
        int m_stations = 0;
        zCTubeSettings m_settings;
        zObjectMesh m_tube;
        zObjectGraph m_sampledDirectrix;
        zObjectMeshArray m_strips, m_unrolled;
    };
}
