#include <zspace/toolsets.h>
#include <iostream>

int main()
{
    // The mesh-walking slicer needs an explicit transverse seed loop. The
    // middle ring has valence four, so its loop stays transverse to the columns.
    zSpace::zObjectMesh mesh;
    zSpace::zFnMesh fn(mesh);
    zSpace::zPointArray positions = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
        {0, 0, 2}, {1, 0, 2}, {1, 1, 2}, {0, 1, 2}
    };
    zSpace::zIntArray counts(10, 4);
    zSpace::zIntArray connects = {
        0, 3, 2, 1, 8, 9, 10, 11,
        0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7,
        4, 5, 9, 8, 5, 6, 10, 9, 6, 7, 11, 10, 7, 4, 8, 11
    };
    fn.create(positions, counts, connects);

    zSpace::zTs3DP threeDP;
    threeDP.fromMesh(mesh);
    threeDP.setPrintLayerHeight(0.25f);
    threeDP.setSlicingSeedVertices(4, 8, 5);
    threeDP.computeAll();

    if (threeDP.sliceMeshes().empty() || threeDP.printPaths().empty())
    {
        std::cerr << "Seeded 3DP fixture produced no slices or print paths.\n";
        return 1;
    }
    for (auto& path : threeDP.printPaths())
    {
        if (zSpace::zFnGraph(path).numEdges() == 0)
        {
            std::cerr << "3DP produced an empty print path.\n";
            return 1;
        }
    }
    return 0;
}
