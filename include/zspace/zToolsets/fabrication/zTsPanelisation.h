#pragma once
#include <zspace/interface.h>
#include <zspace/zToolsets/zToolsetsExport.h>
#include <vector>

namespace zSpace {
struct zPanelisationSettings {
    double spacing = 1.0, offset = 0.0, width = 1.0;
    double seedBand = 0.05, weldTolerance = 0.0001;
    double creaseAngleDegrees = 35.0, smallPanelArea = 0.01;
    int axis = 0; // 0 automatic (longest bounds axis), 1 X, 2 Y, 3 Z
    bool reverse = false, surfaceDistance = true, useCreases = true;
};
struct zPanelisationSource {
    zObjectMesh mesh;
    zIntArray creaseEdges; // vertex endpoint pairs, not edge ordinals
    zIntArray panelSeamSeeds; // mesh edges to extend through opposite edges at regular quad vertices
    bool topologyGuides = false;
    int layoutMode = 0; // 0 existing wall, 1 aligned ceiling, 2 staggered arch
    int alignmentSource = -1;
    std::vector<double> alignmentDirection; // unit horizontal axis for wall stations / arch elevation
    double staggerFraction = 0.5;
    std::vector<double> panelFrame; // optional [Ux,Uy,Uz,Uorigin,Vx,Vy,Vz,Vorigin], in local metres
};
struct zPanelInfo {
    int source = 0, station = 0;
    double area = 0.0;
    int column = 0, chart = 0;
};
struct zPanelisationGuides {
    std::vector<double> positions; // independent segment endpoints, xyz xyz
    zIntArray kinds, sources; // 0 top-seeded blue, 1 inserted red
    int seeds = 0, branches = 0, unsupportedCharts = 0;
    double maxSampledWidth = 0;
};
// Intrinsic quad-row charts with metric X/Y clipping; conformal charts for triangle patches.
// Panel IDs identify connected surface regions, not individual render polygons.
class ZSPACE_TOOLSETS zTsPanelisation {
public:
    void compute(std::vector<zPanelisationSource>& sources, const zPanelisationSettings& settings);
    zObjectMesh& mesh() { return m_mesh; }
    const zIntArray& facePanels() const { return m_facePanels; }
    const zIntArray& faceSources() const { return m_faceSources; }
    const zIntArray& jointEdges() const { return m_jointEdges; }
    const zIntArray& jointKinds() const { return m_jointKinds; } // 0 outer, 1 X, 2 crease/chart/source seam, 3 Y, 4 topology divider
    const std::vector<zPanelInfo>& panels() const { return m_panels; }
    const zPanelisationGuides& guides() const { return m_guides; }
    int components() const { return m_components; }
    int weldedEdges() const { return m_weldedEdges; }
    int smallPanels() const { return m_smallPanels; }
    int resolvedAxis() const { return m_axis; }
private:
    zPanelisationGuides m_guides;
    zObjectMesh m_mesh;
    zIntArray m_facePanels, m_faceSources, m_jointEdges, m_jointKinds;
    std::vector<zPanelInfo> m_panels;
    int m_components = 0, m_weldedEdges = 0, m_smallPanels = 0, m_axis = 0;
};
}
