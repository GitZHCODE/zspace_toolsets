"""Compile the actual Field.cpp edge kernels with lightweight graph stubs.

Run from a Visual Studio developer shell: python scripts/test_field_edges.py
This checks the private kernels, not the full zSpace slicer integration.
"""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if not shutil.which("cl"):
    vcvars = Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)")) / (
        "Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat")
    setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True,
                           check=True, capture_output=True, text=True)
    for line in setup.stdout.splitlines():
        key, sep, value = line.partition("=")
        if sep and key:
            os.environ[key] = value
source = (ROOT / "cpp/source/zToolsets/ts3dp/utils/Field.cpp").read_text()
kernels = source[source.index("struct FastSegment {"):
                 source.index("// Shared by ratio-based and anchored slots")]
prelude = r'''
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
struct zPoint { double x, y, z = 0; };
using zPointArray = std::vector<zPoint>;
using zIntArray = std::vector<int>;
using zScalar = double;
using zScalarArray = std::vector<zScalar>;
struct zObjGraph { zPointArray points; zIntArray edges; };
struct zFnGraph {
  zObjGraph &g;
  explicit zFnGraph(zObjGraph &graph) : g(graph) {}
  void getVertexPositions(zPointArray &p) { p = g.points; }
  void getEdgeData(zIntArray &e) { e = g.edges; }
};
'''
checks = r'''
zObjGraph loop(zPointArray p) {
  zObjGraph g; g.points = p;
  for (int i = 0; i < int(p.size()); ++i) {
    g.edges.push_back(i); g.edges.push_back((i + 1) % p.size());
  }
  return g;
}
// Independent original polygon test, including its endpoint arithmetic.
double referenceSign(const zObjGraph &g, zPoint p) {
  float px = float(p.x), py = float(p.y); int winding = 0;
  for (size_t i = 0; i < g.edges.size(); i += 2) {
    const auto &a = g.points[g.edges[i]], &b = g.points[g.edges[i + 1]];
    float x0 = float(a.x), y0 = float(a.y), x1 = float(b.x), y1 = float(b.y);
    if (py < std::min(y0,y1) || py > std::max(y0,y1)) continue;
    if (y0 <= py) {
      if (y1 > py && (x1-x0)*(py-y0)-(px-x0)*(y1-y0) > 0) ++winding;
    } else if (y1 <= py && (x1-x0)*(py-y0)-(px-x0)*(y1-y0) < 0) --winding;
  }
  return winding ? -1.0 : 1.0;
}
double referenceDistance(const zObjGraph &g, zPoint p) {
  float best = std::numeric_limits<float>::max(); bool valid = false;
  for (size_t i = 0; i < g.edges.size(); i += 2) {
    const auto &a = g.points[g.edges[i]], &b = g.points[g.edges[i+1]];
    float x = float(a.x), y = float(a.y), dx = float(b.x)-x, dy = float(b.y)-y;
    float len = dx*dx+dy*dy; if (len < 1e-12f) continue;
    valid = true;
    float t = ((float(p.x)-x)*dx+(float(p.y)-y)*dy)*(1.0f/len);
    t = std::max(0.0f,std::min(1.0f,t));
    float qx = float(p.x)-(x+t*dx), qy = float(p.y)-(y+t*dy);
    best = std::min(best,qx*qx+qy*qy);
  }
  return valid ? double(std::sqrt(best)) : std::numeric_limits<zScalar>::max();
}
int main() {
  std::vector<zObjGraph> cases = {
    loop({{0,0},{4,0},{4,4},{0,4}}),
    loop({{0,0},{4,0},{4,1},{1,1},{1,4},{0,4}}),
    loop({{0,0},{1e-8,0},{1e-8,1e-8},{0,1e-8}}),
    loop({{0,0},{0,0},{2,0},{2,2},{0,2}}),
    loop({{0,0},{0,0}}), zObjGraph{}
  };
  auto hole = loop({{1,1},{1,3},{3,3},{3,1}});
  auto holed = cases[0];
  for (int v : hole.edges) holed.edges.push_back(v + int(holed.points.size()));
  holed.points.insert(holed.points.end(), hole.points.begin(), hole.points.end());
  cases.push_back(holed);
  auto reversed = cases[1];
  for (size_t i=0;i<reversed.edges.size();i+=2) std::swap(reversed.edges[i],reversed.edges[i+1]);
  cases.push_back(reversed);
  zPointArray points = {{0,0},{4,4},{2,0},{2,2},{1,1},{1e-9,1e-9},{5,5}};
  std::mt19937 rng(42); std::uniform_real_distribution<double> dist(-1,5);
  for (int i=0;i<10000;++i) points.push_back({dist(rng),dist(rng)});
  for (const auto &g : cases) {
    std::vector<FastSegment> edges; computeFastSegments(g,edges);
    zScalarArray signs, distances;
    computeSDF(points.data(),int(points.size()),edges,signs);
    computeDF(points.data(),int(points.size()),edges,distances);
    // Simulate concatenating separately prepared boundary width runs.
    std::vector<FastSegment> partitioned;
    for (size_t e=0;e<g.edges.size();e+=2) {
      zObjGraph run; run.points={g.points[g.edges[e]],g.points[g.edges[e+1]]}; run.edges={0,1};
      std::vector<FastSegment> part; computeFastSegments(run,part);
      partitioned.insert(partitioned.end(),part.begin(),part.end());
    }
    zScalarArray grouped; computeSDF(points.data(),int(points.size()),partitioned,grouped);
    assert(signs==grouped);
    for (size_t i=0;i<points.size();++i) {
      assert(signs[i]==referenceSign(g,points[i]));
      assert(distances[i]==referenceDistance(g,points[i]));
    }
  }
  std::vector<FastSegment> edges; computeFastSegments(holed,edges);
  zPoint probes[]={{0.5,0.5},{2,2},{5,5}}; zScalarArray signs;
  computeSDF(probes,3,edges,signs);
  assert(signs[0]==-1 && signs[1]==1 && signs[2]==1);
  std::cout << "PASS: " << cases.size()*points.size()
            << " distance/sign comparisons; grouped edges and holes\n";
}
'''
with tempfile.TemporaryDirectory(prefix="ts3dp_field_") as tmp:
    cpp = Path(tmp) / "field_edges.cpp"
    exe = Path(tmp) / "field_edges.exe"
    cpp.write_text(prelude + kernels + checks)
    subprocess.run(["cl", "/nologo", "/EHsc", "/std:c++17", "/fp:precise",
                    str(cpp), "/Fe:" + str(exe), "/Fo:" + str(Path(tmp)/"field_edges.obj")],
                   check=True, cwd=tmp)
    subprocess.run([str(exe)], check=True)
