#include <zspace/slicing.h>
#include <zspace/toolsets.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using namespace zSpace;

namespace {
void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

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

void checkPaths(const PrintPathArray& layers) {
  require(layers.size() == 5, "Expected five layers for height 2 and spacing 0.5");
  for (const auto& layer : layers) {
    require(!layer.empty(), "Empty print-path layer");
    for (const auto& path : layer) {
      auto& graph = const_cast<zObjGraph&>(path.graph);
      zFnGraph fn(graph);
      require(fn.numVertices() > 2 && fn.numEdges() > 1, "Degenerate path");
      require(path.widths.size() == static_cast<size_t>(fn.numVertices()), "Missing per-vertex print widths");
      for (float width : path.widths)
        require(std::isfinite(width) && width > 0, "Non-positive or non-finite print width");
      const auto* points = fn.getRawVertexPositions();
      for (int i = 0; i < fn.numVertices(); ++i) {
        const auto& p = points[i];
        require(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z), "Non-finite path position");
        require(p.x >= -0.01 && p.x <= 4.01 && p.y >= -0.01 && p.y <= 4.01 && p.z >= -0.01 && p.z <= 2.01,
                "Path escaped the input block");
      }
    }
  }
}

void topology() {
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  require(session.computeSlicer(false), session.getStatus());
  require(session.getIsoMeshes().size() == 5, "Unexpected slice count");
  require(session.getPathFields().size() == 5, "Missing layer fields");
  require(session.computePrintPath(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  require(session.computePostProcess(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  auto& combined = const_cast<PrintPath&>(session.getResultantPrintPath());
  require(zFnGraph(combined.graph).numEdges() > 0, "PrintSyn did not produce a combined path");

  require(session.computeSlicer(false), session.getStatus());
  require(session.computePostProcess(false), session.getStatus());
  zTsSlicing moved(std::move(session));
  require(moved.computeSlicer(false), moved.getStatus());
  require(moved.computePostProcess(false), moved.getStatus());
  checkPaths(moved.get3DPrintPaths());
  configure(session, mesh);
  require(session.computeSlicer(false), "Moved-from session could not be reused");
}

void invalidInput() {
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  session.clearSlotAnchor();
  require(!session.computeSlicer(false), "Missing anchor was accepted");
  require(session.getStatus().find("Slot") != std::string::npos, "Missing anchor diagnostic");
  session.setSlotAnchor(zPoint(0,0,0));
  session.setInputMesh(mesh, zIntArray(1, Bottom));
  require(!session.computeSlicer(false), "Invalid face tags were accepted");
  require(session.get3DPrintPaths().empty(), "Failure retained stale print paths");
  zTsIO::SlicerUSDInput usdInput;
  usdInput.Block = mesh;
  std::string ioError;
  require(!zTsIO::readSlicerInputUSD("zspace-missing-slicer-input.usda", usdInput, &ioError), "Missing USD input was accepted");
  require(!ioError.empty() && zFnMesh(usdInput.Block).numVertices() == 0, "USD failure did not clear data and report a diagnostic");
  require(!zTsSlicing::supportsSlicerMode(static_cast<zTsSlicerMode>(999)), "Unknown mode reported supported");
  bool rejected = false;
  try { session.setSlicerMode(static_cast<zTsSlicerMode>(999)); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected && session.getSlicerMode() == zTsSlicerMode::Topology, "Mode validation changed the session");
#ifndef ZSPACE_TOOLSETS_WITH_TETGEN
  rejected = false;
  try { session.setSlicerMode(zTsSlicerMode::Harmonic); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "Disabled harmonic mode was accepted");
#endif
}

void features() {
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  zIntArray tags(10, Side), widthIds(10, 0);
  tags[0] = Bottom; tags[1] = Top;
  widthIds[2] = widthIds[6] = 1;
  session.setInputMesh(mesh, tags, widthIds);
  session.getSettings().BoundaryWidth = {0.4f, 0.6f};
  session.getSettings().DoubleBoundary = true;
  session.getSettings().FieldCellSize = 0.05f;
  zObjGraph bracing;
  zPointArray points = {{0,2,0},{4,2,0}};
  zIntArray edges = {0,1};
  zFnGraph(bracing).create(points, edges);
  session.setBracingGraph(bracing);
  require(session.getWidthGroupCount() == 2, "Width groups were lost");
  require(session.computeSlicer(false), session.getStatus());
  require(session.getTopologySlicer()->getBracingGraphs2D().size() == 5, "Bracing was not mapped to each layer");
  const auto& widths = session.getWidthFields();
  bool varied = false;
  for (const auto& layer : widths) {
    if (layer.empty()) continue;
    for (double width : layer) if (std::abs(width - layer.front()) > 1e-3) varied = true;
  }
  require(varied, "Per-group print widths were lost in field computation");
  require(session.computePrintPath(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
}

void parallelParity() {
  auto mesh = block();
  zTsSlicing serial, parallel;
  configure(serial, mesh);
  configure(parallel, mesh);
  require(serial.computeSlicer(false), serial.getStatus());
  require(parallel.computeSlicer(true), parallel.getStatus());
  require(serial.getPathFields() == parallel.getPathFields(), "Serial and parallel path fields differ");
  require(serial.getWidthFields() == parallel.getWidthFields(), "Serial and parallel width fields differ");
  require(serial.computePostProcess(false), serial.getStatus());
  require(parallel.computePostProcess(true), parallel.getStatus());
  const auto& a = serial.get3DPrintPaths();
  const auto& b = parallel.get3DPrintPaths();
  require(a.size() == b.size(), "Serial and parallel layer counts differ");
  for (size_t layer = 0; layer < a.size(); ++layer) {
    require(a[layer].size() == b[layer].size(), "Serial and parallel path counts differ");
    for (size_t i = 0; i < a[layer].size(); ++i) {
      zPointArray pa, pb;
      zFnGraph(const_cast<zObjGraph&>(a[layer][i].graph)).getVertexPositions(pa);
      zFnGraph(const_cast<zObjGraph&>(b[layer][i].graph)).getVertexPositions(pb);
      require(pa.size() == pb.size(), "Serial and parallel vertex counts differ");
      require(a[layer][i].widths == b[layer][i].widths, "Serial and parallel widths differ");
      for (size_t j = 0; j < pa.size(); ++j)
        require(zTsUtils::vecDistSq(pa[j], pb[j]) < 1e-12, "Serial and parallel positions differ");
    }
  }
}

void harmonic() {
#ifdef ZSPACE_TOOLSETS_WITH_TETGEN
  auto mesh = block();
  zTsSlicing session;
  configure(session, mesh);
  session.setSlicerMode(zTsSlicerMode::Harmonic);
  require(session.getHarmonicSlicer() != nullptr, "Harmonic operator unavailable");
  require(session.computeSlicer(false), session.getStatus());
  require(session.computePrintPath(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
  require(session.computePostProcess(false), session.getStatus());
  checkPaths(session.get3DPrintPaths());
#endif
}
}

int main(int argc, char** argv) {
  try {
    const std::string test = argc > 1 ? argv[1] : "topology";
    if (test == "topology") topology();
    else if (test == "invalid") invalidInput();
    else if (test == "features") features();
    else if (test == "parallel") parallelParity();
    else if (test == "harmonic") harmonic();
    else throw std::runtime_error("Unknown test case");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
