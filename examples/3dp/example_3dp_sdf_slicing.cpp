#include <zspace/interface.h>
#include <zspace/io.h>
#include <zspace/toolsets.h>

int main()
{
    zSpace::zTs3DP ThreeDP;
    ThreeDP.setPrintLayerHeight(0.01f);
    ThreeDP.setFieldResolution(200, 80);
    ThreeDP.setPrintWidth(0.004f);
    ThreeDP.setPrintSpacing(0.005f);

    if (!ThreeDP.readMesh("data/3dp/input_quad_mesh.obj")) return 1;

    ThreeDP.computeSlices();
    ThreeDP.computeUnrolledSDFs();
    ThreeDP.computePrintPaths();
    ThreeDP.computePrintMesh();

    return 0;
}
