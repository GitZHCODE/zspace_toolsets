#include <zspace/slicing.h>
#include <iostream>

using namespace zSpace;

namespace {
zObjectMesh block() {
  zObjectMesh mesh;
  zPointArray points = {
    {0,0,0},{4,0,0},{4,4,0},{0,4,0},
    {0,0,1},{4,0,1},{4,4,1},{0,4,1},
    {0,0,2},{4,0,2},{4,4,2},{0,4,2}
  };
  zIntArray counts(10, 4);
  zIntArray vertices = {
    0,3,2,1, 8,9,10,11,
    0,1,5,4, 1,2,6,5, 2,3,7,6, 3,0,4,7,
    4,5,9,8, 5,6,10,9, 6,7,11,10, 7,4,8,11
  };
  zFnMesh(mesh).create(points, counts, vertices);
  return mesh;
}

void configure(zTsSlicing& session, zObjectMesh& mesh) {
  zIntArray tags(10, Side);
  tags[0] = Bottom;
  tags[1] = Top;
  session.setInputMesh(mesh, tags);
  session.setSlotAnchor(zPoint(0, 0, 0));
  auto& settings = session.getSettings();
  settings.LayerSpacing = 0.5f;
  settings.BoundaryWidth = {0.4f};
  settings.FieldCellSize = 0.1f;
  session.getPrintSynSetting().ResampleLength = 0.2f;
}

}

int main() {
  auto mesh = block();
  zTsSlicing slicer;
  configure(slicer, mesh);
  if (!slicer.computeSlicer(false) || !slicer.computePostProcess(false)) {
    std::cerr << slicer.getStatus() << '\n';
    return 1;
  }
  const auto& path = slicer.getResultantPrintPath();
  std::cout << "Layers: " << slicer.get3DPrintPaths().size()
            << "; print length: " << slicer.getPrintLength() << '\n';
  return zFnGraph(const_cast<zObjGraph&>(path.graph)).numEdges() > 0 ? 0 : 1;
}
