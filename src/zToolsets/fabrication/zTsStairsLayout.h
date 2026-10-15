#pragma once
#include <zspace/zToolsets/fabrication/zTsParametricStairs.h>
#include <array>

namespace zSpace {
struct zStairsFootprint { zPointArray polygon; double underside; double top=0; int flight=-1; };
struct zStairsLayoutData {
    zObjectMesh mesh;
    zObjectGraph graph;
    zPointArray graphPoints;
    zIntArray faceTypes, riserCounts, edgeFlights, meshFaceFlights;
    std::vector<std::array<double,2>> flightFloors;
    zDoubleArray actualRisers, usableRuns, startReserves, endReserves;
    std::vector<zIntArray> flightPaths;
    std::vector<zStairsFootprint> footprints;
    double minimumClearance = 0;
};
ZSPACE_TOOLSETS std::shared_ptr<zStairsLayoutData> buildStairsLayout(const zPointArray& points,
    const zIntArray& edges, const zParametricStairsSettings& settings,
    const zIntArray& lockedCounts = {}, bool geometry = true, bool capacityOnly = false);
// Exact projected triangle/polygon overlap extrema, including edge crossings.
ZSPACE_TOOLSETS double stairsClearance(const zPointArray& points, const zIntArray& counts,
    const zIntArray& connects, const zStairsLayoutData& layout, double halfThickness);
ZSPACE_TOOLSETS void lowerStairsBelowEnvelope(zPointArray& points, const zIntArray& counts,
    const zIntArray& connects, const zIntArray& supports,
    const zStairsLayoutData& layout, double halfThickness, double clearance);
}
