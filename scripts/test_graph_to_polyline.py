"""Compile and check the actual graphtopolyline implementation with graph stubs.

Run with py -3 scripts/test_graph_to_polyline.py (MSVC required).
This does not build the full zSpace integration.
"""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
if not shutil.which("cl"):
    vcvars = Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)")) / (
        "Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat")
    setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True,
                           check=True, capture_output=True, text=True)
    for line in setup.stdout.splitlines():
        key, sep, value = line.partition("=")
        if sep and key:
            os.environ[key] = value
source = (root / "cpp/source/zToolsets/ts3dp/utils/Utils.cpp").read_text()
function = source[source.index("ZSPACE_INLINE bool graphtopolyline("):
                  source.index("ZSPACE_INLINE bool samplePolyline(const zObjGraph")]
field_source = (root / "cpp/source/zToolsets/ts3dp/utils/Field.cpp").read_text()
slot = field_source[field_source.index("ZSPACE_INLINE bool createSlot("):
                    field_source.index("ZSPACE_INLINE float widthOfGroup(")]
slot_helper = field_source[field_source.index("// Shared by ratio-based and anchored slots"):
                           field_source.index("} // anonymous namespace")]
prelude = r'''
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#define ZSPACE_INLINE inline
struct zPoint {
  double x,y,z;
  zPoint(double a=0,double b=0,double c=0):x(a),y(b),z(c){}
};
struct zVector {
  double x=0,y=0,z=0;
  double length() const{return std::hypot(x,y,z);}
  void normalize(){double l=length(); x/=l;y/=l;z/=l;}
};
zVector vecSub(const zPoint&a,const zPoint&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
using zPointArray=std::vector<zPoint>;
using zIntArray=std::vector<int>;
using zScalar=double;
using zScalarArray=std::vector<zScalar>;
struct zObjGraph { zPointArray points; zIntArray edges; };
struct zFnGraph {
  zObjGraph &g;
  explicit zFnGraph(zObjGraph &graph):g(graph){}
  void getVertexPositions(zPointArray &p){p=g.points;}
  void getEdgeData(zIntArray &e){e=g.edges;}
  void clear(){g.points.clear();g.edges.clear();}
  void create(zPointArray &p,zIntArray &e){g.points=p;g.edges=e;}
};
bool samplePolyline(const zPoint*,const zScalar*,size_t,double,zPoint&,
                    zVector*,double minSegmentLength=1e-6);
'''
checks = r'''
void checkSlot(const zObjGraph &g,double t,zPoint center,double halfLength) {
  zObjGraph result; assert(createSlot(g,result,t,halfLength));
  assert(result.points.size()==2 && result.edges==zIntArray({0,1}));
  const auto &a=result.points[0], &b=result.points[1];
  assert(std::abs((a.x+b.x)*0.5-center.x)<1e-10);
  assert(std::abs((a.y+b.y)*0.5-center.y)<1e-10);
  assert(a.z==center.z && b.z==center.z);
  assert(std::abs(std::hypot(a.x-b.x,a.y-b.y)-2*halfLength)<1e-10);
}
int main() {
  zPointArray p; zScalarArray arc;
  zObjGraph open{{{2,0},{3,0},{0,0},{1,0}}, {0,3,1,0,3,2}};
  assert(graphtopolyline(open,p,arc));
  assert(p.size()==4 && arc==zScalarArray({0,1,2,3}));
  for(int i=0;i<4;++i) assert(p[i].x==3-i);
  assert(open.points[0].x==2 && open.edges==zIntArray({0,3,1,0,3,2}));
  zObjGraph closed{{{0,0},{2,0},{2,2},{0,2}}, {2,3,0,1,3,0,1,2}};
  assert(graphtopolyline(closed,p,arc));
  assert(p.size()==5 && arc==zScalarArray({0,2,4,6,8}));
  assert(p.front().x==p.back().x && p.front().y==p.back().y);
  zObjGraph line3D{{{0,0,0},{0,0,3},{0,4,3}}, {1,2,0,1}};
  assert(graphtopolyline(line3D,p,arc));
  assert(arc==zScalarArray({0,3,7}));
  zObjGraph tiny{{{0,0},{1e-10,0}}, {0,1}};
  assert(graphtopolyline(tiny,p,arc) && arc.back()==1e-10);
  const double nan=std::numeric_limits<double>::quiet_NaN();
  const double inf=std::numeric_limits<double>::infinity();
  std::vector<zObjGraph> invalid={
    {}, {{{0,0}},{}}, {{{0,0},{1,0}},{}},
    {{{0,0},{1,0}},{0}},
    {{{0,0},{1,0}},{0,2}}, {{{0,0},{1,0}},{-1,1}},
    {{{0,0},{1,0}},{0,0}}, {{{0,0},{1,0}},{0,1,1,0}},
    {{{0,0},{0,0}},{0,1}},
    {{{0,0},{1,0},{2,0}},{0,1}},
    {{{0,0},{1,0},{0,1},{-1,0}},{0,1,0,2,0,3}},
    {{{0,0},{1,0},{2,0},{3,0}},{0,1,2,3}},
    {{{0,0},{1,0},{0,1},{3,0},{4,0},{3,1}},{0,1,1,2,2,0,3,4,4,5,5,3}},
    // An open chain plus a separate loop passes the edge-count check.
    {{{0,0},{1,0},{3,0},{4,0},{3,1}},{0,1,2,3,3,4,4,2}},
    {{{nan,0},{1,0}},{0,1}}, {{{0,0},{inf,0}},{0,1}},
    {{{-1e308,0},{1e308,0}},{0,1}}
  };
  for(const auto &g:invalid) {
    p={{9,9}}; arc={9};
    assert(!graphtopolyline(g,p,arc));
    assert(p.empty() && arc.empty());
    zObjGraph output{{{9,9},{10,10}},{0,1}};
    assert(!createSlot(g,output,0.5,1));
    assert(output.points.empty() && output.edges.empty());
  }
  checkSlot(open,0.5,{1.5,0},2);
  checkSlot(open,-1,{3,0},1);
  checkSlot(open,2,{0,0},1);
  checkSlot(closed,0.25,{2,0},1);
  checkSlot(closed,1,{0,0},1);
  // 3D arc length sampling, but a full-length XY perpendicular.
  zObjGraph sloped{{{0,0,0},{3,0,4}},{0,1}};
  checkSlot(sloped,0.5,{1.5,0,2},2);
  zObjGraph output;
  assert(!createSlot(line3D,output,0,1)); // vertical tangent
  assert(!createSlot(open,output,nan,1));
  assert(!createSlot(open,output,0.5,inf));
  assert(!createSlot(open,output,0.5,0));
  assert(!createSlot(open,output,0.5,-1));
  auto alias=open;
  assert(createSlot(alias,alias,0.5,1) && alias.points.size()==2);
  std::cout<<"PASS: 4 valid polylines and "<<invalid.size()
           <<" rejected graphs; slot centers, half-length, clamps, invalid inputs and aliasing\n";
}
'''
with tempfile.TemporaryDirectory(prefix="ts3dp_polyline_") as tmp:
    cpp = Path(tmp) / "graph_to_polyline.cpp"
    exe = Path(tmp) / "graph_to_polyline.exe"
    cpp.write_text(prelude + function + slot_helper + slot + checks)
    subprocess.run(["cl", "/nologo", "/EHsc", "/std:c++17", str(cpp),
                    "/Fe:" + str(exe), "/Fo:" + str(Path(tmp)/"graph_to_polyline.obj")],
                   check=True, cwd=tmp)
    subprocess.run([str(exe)], check=True)
