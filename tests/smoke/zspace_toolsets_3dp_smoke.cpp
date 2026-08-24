#include <zspace/toolsets.h>

int main()
{
    zSpace::zObjectMesh mesh;
    zSpace::zFnMesh fn(mesh);

    zSpace::zPointArray positions = {
        zSpace::zPoint(0, 0, 0),
        zSpace::zPoint(1, 0, 0),
        zSpace::zPoint(1, 1, 0),
        zSpace::zPoint(0, 1, 0),
        zSpace::zPoint(0, 0, 1),
        zSpace::zPoint(1, 0, 1),
        zSpace::zPoint(1, 1, 1),
        zSpace::zPoint(0, 1, 1)
    };
    zSpace::zIntArray counts = { 4, 4, 4, 4, 4, 4 };
    zSpace::zIntArray connects = {
        0, 1, 2, 3,
        4, 7, 6, 5,
        0, 4, 5, 1,
        1, 5, 6, 2,
        2, 6, 7, 3,
        3, 7, 4, 0
    };
    fn.create(positions, counts, connects);

    zSpace::zTs3DP ThreeDP;
    ThreeDP.fromMesh(mesh);
    ThreeDP.setPrintLayerHeight(0.25f);
    ThreeDP.computeAll();

    return ThreeDP.printPaths().empty() ? 1 : 0;
}
