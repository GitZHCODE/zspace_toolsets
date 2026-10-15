#pragma once
#include <zspace/interface.h>
#include <zspace/zToolsets/zToolsetsExport.h>
#include <memory>
#include <string>
#include <functional>

namespace zSpace {
enum class zVaultPlanConstraint { Free = 0, PenalizeXY = 1, FixXY = 2 };
struct zVaultDiagramSettings { double tolerance = 1e-5; };
struct zVaultDiagramReport {
    std::string status = "uninitialized";
    double horizontalResidual = 0, closureError = 0, maxEdgeMismatch = 0;
    bool exact = false;
};
struct zVaultBestFitSettings {
    // Magnitudes: tension [qMin,qMax], compression [-qMax,-qMin].
    double qMin = 0.01, qMax = 1000.0;
    bool compressionOnly = true;
    bool fixXY = false; // Preserve input X/Y; only free-vertex Z may move.
    zVaultPlanConstraint planConstraint = zVaultPlanConstraint::Free;
    double xyPenaltyWeight = 1.0;
    double geometryWeight = 0.5;
    double minMaxEdgeRatio = 0.1; // Desired minimum / maximum active 3D edge length, [0,1].
    double edgeRegularizationStrength = 0.0; // Zero disables regularization.
    double absoluteTolerance = 1e-5, relativeTolerance = 0.01;
    double initialStep = 0.1, integrationTolerance = 1e-4;
    int maxIterations = 2000;
};
struct zVaultSolveReport {
    std::string status = "uninitialized";
    int iteration = 0, rank = 0, activeEdges = 0, freeVertices = 0; // rank=-1: sparse iterative fit, not measured.
    double residualNorm = 0, maxResidual = 0, initialResidual = 0;
    double maxDisplacement = 0, stepSize = 0;
    double minMaxEdgeRatio = 0, edgeRegularizationError = 0;
    double rmsXYDisplacement = 0, maxXYDisplacement = 0;
    bool converged = false;
};
// Lachauer & Block (2014), Interactive Equilibrium Modelling, sections 3.1-3.4.
// Geometry and assignments are snapshots. Setters invalidate the current solve.
// q > 0 is tension. Loads = explicit loads + mass * gravity, in consistent units.
class ZSPACE_TOOLSETS zTsVault {
public:
    zTsVault();
    ~zTsVault();
    zTsVault(const zTsVault&) = delete;
    zTsVault& operator=(const zTsVault&) = delete;
    void setInputMesh(zObjectMesh& mesh);
    // Separate attraction reference; retained through settings, fit, reset and rewind.
    void setReferencePositions(const zPointArray& positions);
    // Candidate feasibility hook: rejecting a trial leaves current geometry intact.
    void setGeometryValidator(std::function<bool(const zPointArray&)> validator);
    void setSupports(const zIntArray& vertexIds);
    void setTensionEdges(const zIntArray& edgeIds);
    void setVertexMass(double mass);
    void setVertexMasses(const zDoubleArray& masses);
    void setGravity(const zVector& gravity);
    void setVertexLoads(const zVectorArray& loads);
    // Optional signed bounds, one pair per original mesh edge; sign assignments
    // are intersected with these bounds. Empty arrays restore settings bounds.
    void setForceDensityBounds(const zDoubleArray& lower, const zDoubleArray& upper);
    void setBestFitSettings(const zVaultBestFitSettings& settings);
    zVaultSolveReport computeBestFitForceDensities();
    zVaultSolveReport stepBestFit(int iterations = 1);
    zVaultSolveReport solveBestFit(int maxIterations);
    // Direct equilibrium for current bounded densities, respecting fixed XY
    // when enabled; full XYZ residuals remain the convergence criterion.
    zVaultSolveReport snapFixedDensities(double geometryWeight = 0.0);
    void resetBestFit();
    bool previousBestFit(); // Restores one accepted step, including densities.
    void clear();
    zObjectMesh& resultMesh(); // Treat as read-only while solving.
    const zDoubleArray& forceDensities() const;
    const zDoubleArray& axialForces() const;
    const zVectorArray& residualForces() const;
    const zVectorArray& supportReactions() const;
    const zIntArray& activeEdgeIds() const;
    const zIntArray& edgeConnects() const;
    const zVaultSolveReport& report() const;
    // Stage 1: reconstruct diagrams without changing geometry or force densities.
    zVaultDiagramReport createFormFromResult();
    zVaultDiagramReport createForceFromForm(const zVaultDiagramSettings& settings = {});
    zVaultDiagramReport updateForceFromDensities();
    const zObjectMesh& formMesh() const;
    const zObjectMesh& forceMesh() const;
    const zIntArray& formHalfEdgeToForceHalfEdge() const;
    const zIntArray& forceHalfEdgeToFormHalfEdge() const;
    const zVaultDiagramReport& diagramReport() const;
private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};
}
