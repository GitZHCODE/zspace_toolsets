#pragma once
#include <zspace/interface.h>
#include <zspace/zToolsets/zToolsetsExport.h>
#include <memory>
#include <string>

namespace zSpace {
enum class zForcePolyhedraFaceType { Member, Load, Reaction, Zero, Virtual };
struct zForcePolyhedraSource {
    zForcePolyhedraFaceType type = zForcePolyhedraFaceType::Member;
    int id = -1; // Original graph edge for Member; vertex for Load/Reaction.
    double magnitude = 0;
};
struct zForcePolyhedraFace {
    zVector normal;
    double targetArea = 0, area = 0;
    zForcePolyhedraFaceType type = zForcePolyhedraFaceType::Zero;
    // Parallel, same-direction forces share a face; no source is discarded.
    std::vector<zForcePolyhedraSource> sources;
};
struct zForcePolyhedraSettings {
    double equilibriumTolerance = 1e-6;
    double areaTolerance = 1e-7; // Relative to each physical target; zero faces use the largest target.
    double correspondingAreaTolerance = 1e-4; // Absolute |A0-A1| for corresponding member-face contributions.
    double normalTolerance = 1e-7; // Sine of normal-angle error; also bounds planarity / cell size.
    double forceScale = 1.0; // Area = force magnitude * forceScale.
    int maxIterations = 5000;
    // Movement weight: 0 = force cells only; 1 = form graph only.
    double formWeight = 0.5;
    double angleTolerance = 0.01; // Degrees: edge/face and opposite contact normals.
    double relaxation = 0.5;
    double minMaxEdgeRatio = 0.01; // Legacy Polytopal length clamp.
    int areaIterationsPerStep = 10;
    double timeStep = 1.0; // Overdamped unit-mass Euler integration.
    double maxVertexStepRatio = 0.05; // Maximum step / initial diagram extent.
    double normalStepDegrees = 5.0; // Maximum face-normal goal rotation per pass.
};
struct zForcePolyhedraReport {
    std::string status = "uninitialized";
    bool equilibrium = false, cellsConverged = false;
    bool areaPursuitResolved = false; // Corresponding areas pass independently of normal/planarity errors.
    double maxCorrespondingAreaDifference = 0;
    int worstAreaDifferenceEdge = -1;
    int iterations = 0;
    double maxResidual = 0, residualNorm = 0;
    double maxAreaError = 0, maxNormalError = 0, maxPlanarityError = 0;
    double maxAngleDeviation = 0, maxContactAreaMismatch = 0;
    int equilibriumIterations = 0;
    double maxCellVertexStep = 0, maxFormVertexStep = 0;
    int worstAreaCell = -1, worstAngleEdge = -1, worstMotionCell = -1;
};
struct zForcePolyhedraCellReport {
    std::string status = "uninitialized";
    bool converged = false;
    int forceRank = 0, iterations = 0;
    double closureError = 0, maxAreaError = 0, maxNormalError = 0, maxPlanarityError = 0;
    double maxVertexStep = 0;
    int dynamicsSteps = 0;
};

// Lee, Van Mele & Block (2018), Disjointed force polyhedra, sections 3.2-3.4.
// Graph snapshot + one independent mesh per node. Inputs are in consistent units.
// Positive axial force means tension. Loads include mass * gravity.
// EGI and AP do not imply structural stability or constitutive/material analysis.
class ZSPACE_TOOLSETS zTsForcePolyhedra {
public:
    zTsForcePolyhedra();
    ~zTsForcePolyhedra();
    zTsForcePolyhedra(const zTsForcePolyhedra&) = delete;
    zTsForcePolyhedra& operator=(const zTsForcePolyhedra&) = delete;
    void clear();
    void setInputGraph(zObjectGraph& graph); // Resets assignments, settings and results.
    // Edit the form after AP without reconstructing the force diagram. Pinned
    // vertices must retain their input coordinates. Enables force-fixed solving.
    void setFormPositions(const zPointArray& positions);
    void setSupports(const zIntArray& vertexIds); // Fully pinned XYZ supports.
    void setTensionEdges(const zIntArray& edgeIds); // Other members are compression.
    // Positive target force magnitudes, one per graph edge. Empty resets to an
    // automatic seed based on the nodal load scale. AP preserves these targets.
    void setEdgeForceMagnitudes(const zDoubleArray& magnitudes);
    void setVertexMass(double mass);
    void setVertexMasses(const zDoubleArray& masses);
    void setGravity(const zVector& gravity);
    void setVertexLoads(const zVectorArray& loads); // Added to self weight, including at supports.
    void setSettings(const zForcePolyhedraSettings& settings);
    void setFormWeight(double weight); // Continue the current session; no reconstruction.
    // Build EGI arc arrangement + dual seed meshes, then adaptive-normal AP.
    // Planar/collinear nodes receive untargeted virtual
    // faces (paper Fig. 6); these never enter equilibrium or physical-force totals.
    zForcePolyhedraReport createForcePolyhedra();
    zForcePolyhedraReport stepAreaPursuit(int iterations = 1);
    // 3DGS reciprocation, reimplemented from legacy zTsPolytopal geometric
    // targets and neighbour-projection relaxation. Couples each graph member
    // to BOTH corresponding cell-face normals; preserves supports and loads.
    zForcePolyhedraReport stepEquilibrium(int iterations = 1);
    zForcePolyhedraReport solveEquilibrium(int maxIterations = 2000);
    void reset(); // Restore original graph, retain loads/supports/settings.
    const zObjectGraph& resultGraph() const;
    const zIntArray& edgeConnects() const;
    const zDoubleArray& axialForces() const;
    const zVectorArray& residualForces() const;
    const zVectorArray& supportReactions() const;
    // Cell geometry is centred at the origin, one cell per original node.
    // Empty mesh with status=no_forces for unloaded zero-force nodes.
    const zObjectMesh& forceCell(int vertexId) const;
    const std::vector<zForcePolyhedraFace>& cellFaces(int vertexId) const;
    const zForcePolyhedraCellReport& cellReport(int vertexId) const;
    const zForcePolyhedraReport& report() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
}
