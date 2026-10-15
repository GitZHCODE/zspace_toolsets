#pragma once
#include <zspace/zToolsets/statics/zTsVault.h>

namespace zSpace {
struct zStairsLayoutData;
struct zParametricStairsSettings {
    double width = 1.0; // Full ribbon width, in model units.
    int numDivs = 1; // Face subdivision levels; zero keeps the coarse mesh.
    int numSmoothDivs = 1; // Catmull-Clark levels BEFORE Best Fit FDM.
    int divisions = -1; // Legacy longitudinal sampling, when explicitly >= 2.
    double miterLimit = 8.0; // Reject sharp junctions beyond this half-width multiple.
    double vertexMass = 1.0;
    bool generateSteps = false; // Legacy ribbon-only inputs remain compatible.
    double treadDepth = .30, riserHeight = .17, landingDepth = 1.0;
    double stepThickness = .08, structureThickness = .15, clearance = .05;
    zIntArray landingVertices; // Stable landing-node classification for selected layouts.
    zIntArray flightRiserCounts; // Optional per-flight layout selected by graph optimisation.
    zDoubleArray flightTreadDepths;
    std::vector<zPointArray> landingPolygons; // Horizontal closed Rhino outlines.
    zDoubleArray vertexWidths; // Optional full widths, indexed by graph vertex.
    zVaultBestFitSettings bestFit = [] {
        zVaultBestFitSettings s; s.compressionOnly=true;
        s.fixXY=true; s.planConstraint=zVaultPlanConstraint::FixXY; return s;
    }();
};
struct zStairsGraphOptimisation {
    zObjectGraph graph;
    std::vector<zPointArray> landingPolygons;
    zDoubleArray landingOffsets, treadDepths, riserHeights;
    zIntArray riserCounts, bottomNodes, landingVertices;
    int adjustedLandings = 0;
};
// Graph is an XY-plan ribbon layout with arbitrary elevations. Geometry setters
// regenerate a snapshot and initialise Best Fit FDM; stepBestFit advances it.
// End-cap centre and side vertices inherit supports from degree-one graph nodes.
class ZSPACE_TOOLSETS zTsParametricStairs {
public:
    static zStairsGraphOptimisation optimiseGraph(zObjectGraph& graph,
        const zParametricStairsSettings& settings, double treadMin, double treadMax,
        double riserMin, double riserMax);
    void update(zObjectGraph& graph, const zParametricStairsSettings& settings = {});
    zVaultSolveReport stepBestFit(int iterations = 1);
    zVaultSolveReport snapFixedDensities(double geometryWeight = 0.0);
    void smoothResult(int divisions = 1);
    void resetBestFit();
    bool previousBestFit();
    zObjectMesh& centerMesh() { return m_mesh; }
    zObjectMesh& lowPolyMesh() { return m_lowPoly; }
    zObjectMesh& subdividedMesh() { return m_subdivided; }
    zObjectMesh& resultMesh();
    bool hasStairLayout() const { return bool(m_layout); }
    zObjectMesh& stairMesh();
    zObjectGraph& walkingGraph();
    const zIntArray& stairFaceTypes() const;
    const zIntArray& flightRiserCounts() const;
    const zDoubleArray& flightActualRisers() const;
    double minimumClearance() const;
    void refreshLayout();
    const zIntArray& supports() const { return m_supports; }
    const zIntArray& meshVertexGraphVertices() const { return m_sources; }
    const zIntArray& centerGraphMeshVertices() const { return m_centers; }
    const zDoubleArray& forceDensities() const;
    const zIntArray& edgeConnects() const;
    const zIntArray& activeEdgeIds() const;
    const zDoubleArray& axialForces() const;
    const zVectorArray& residualForces() const;
    const zVectorArray& supportReactions() const;
    const zVaultSolveReport& report() const;
private:
    zObjectMesh m_mesh, m_lowPoly, m_subdivided;
    zIntArray m_supports, m_sources, m_centers;
    std::unique_ptr<zTsVault> m_vault;
    zParametricStairsSettings m_settings;
    std::shared_ptr<zStairsLayoutData> m_layout;
    zPointArray m_designPoints, m_centerReferences;
    zIntArray m_designEdges, m_flightCounts;
    std::shared_ptr<zStairsLayoutData> layoutFor(const zPointArray& meshPoints, bool geometry = false) const;
    void attachClearanceValidator();
};
}
