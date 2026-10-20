"""Compare actual field kernels against a pre-refactor source snapshot.

Usage: py -3 scripts/test_field_composition.py path/to/old/Field.cpp
Uses graph/grid stubs and max(a,-b) for boolean_subtract(false); this
checks composition equivalence, not the unavailable zSpace core binary.
"""
import ast
import os
from pathlib import Path
import shutil
import subprocess
import sys
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
old = Path(sys.argv[1]).read_text()
new = (root / "cpp/source/zToolsets/ts3dp/utils/Field.cpp").read_text()
utils = (root / "cpp/source/zToolsets/ts3dp/utils/Utils.cpp").read_text()
tree = ast.parse((root / "scripts/test_graph_to_polyline.py").read_text())
prelude = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign)
               and any(isinstance(t, ast.Name) and t.id == "prelude" for t in n.targets))
prelude = prelude.replace("void clear(){", "int numEdges(){return int(g.edges.size()/2);}\n  void clear(){")
stubs = r'''
#include <random>
struct zObjMeshScalarField { zPointArray points; };
struct zFnMesh {
  zObjMeshScalarField &f;
  explicit zFnMesh(zObjMeshScalarField &field):f(field){}
  int numVertices(){return int(f.points.size());}
  const zPoint *getRawVertexPositions(){return f.points.empty()?nullptr:f.points.data();}
};
struct zFnMeshScalarField {
  zObjMeshScalarField &f;
  explicit zFnMeshScalarField(zObjMeshScalarField &field):f(field){}
  int numFieldValues(){return int(f.points.size());}
  void boolean_subtract(const zScalarArray&a,const zScalarArray&b,zScalarArray&out,bool){
    assert(a.size()==b.size());out.resize(a.size());
    for(size_t i=0;i<a.size();++i)out[i]=std::max(a[i],-b[i]);
  }
};
struct FieldSegment {zObjGraph graph;int type=0;zScalar width=0,offset=0,trim=0;};
struct LoopAnchor {
  int edge=-1;double t=0;zPoint point;bool ccw=true;double distance=0;
  bool valid()const{return edge>=0;}
};
double vecDistSq(const zPoint&a,const zPoint&b){
  return (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z);
}
void combineGraphObjects(const std::vector<zObjGraph>&graphs,zObjGraph&out){
  out={};for(const auto&g:graphs){int base=int(out.points.size());
    out.points.insert(out.points.end(),g.points.begin(),g.points.end());
    for(int v:g.edges)out.edges.push_back(v+base);
  }
}
'''
def implementations(s):
    return (s[s.index("ZSPACE_INLINE bool createSlot("):s.index("ZSPACE_INLINE int loopWidthSegments(")]
            + s[s.index("ZSPACE_INLINE bool createAnchoredSlot("):s.index("ZSPACE_INLINE bool locateSlotAnchor(")]
            + s[s.index("ZSPACE_INLINE void\ncomputeSegmentDistances("):s.index("ZSPACE_INLINE void\nextractRawIsoContour(")])
shared = (utils[utils.index("ZSPACE_INLINE bool graphtopolyline("):
                utils.index("ZSPACE_INLINE bool samplePolyline(const zObjGraph")]
          + new[new.index("struct FastSegment {"):new.index("} // anonymous namespace")])
checks = r'''
void equal(const zScalarArray&a,const zScalarArray&b){
  assert(a.size()==b.size());for(size_t i=0;i<a.size();++i)
    assert(std::abs(a[i]-b[i])<=1e-10*(1+std::abs(a[i])));
}
void equalGraph(const zObjGraph&a,const zObjGraph&b){
  assert(a.edges==b.edges && a.points.size()==b.points.size());
  for(size_t i=0;i<a.points.size();++i)assert(vecDistSq(a.points[i],b.points[i])<1e-20);
}
int main(){
  zObjMeshScalarField grid;
  std::mt19937 rng(23);std::uniform_real_distribution<double>d(-1,5);
  for(int i=0;i<2000;++i)grid.points.push_back({d(rng),d(rng)});
  zObjGraph boundary{{{0,0},{4,0},{4,4},{0,4}},{0,1,1,2,2,3,3,0}};
  std::vector<FieldSegment> base;
  for(int e=0;e<4;++e){FieldSegment s;
    s.graph={{boundary.points[e],boundary.points[(e+1)%4]},{0,1}};
    s.width=0.5+0.2*e;s.offset=0.5*s.width;base.push_back(s);
  }
  FieldSegment rib;rib.type=1;rib.width=0.6;rib.trim=0.4;
  rib.graph={{{0.5,0.7},{2,0.7},{3.5,1.2}},{0,1,1,2}};
  auto rib2=rib;rib2.trim=0;rib2.width=0.8;
  rib2.graph={{{1,0.5},{1,3.5}},{0,1}};
  zObjGraph slot{{{2,-1},{2,1}},{0,1}};
  int comparisons=0;
  for(int mode=0;mode<4;++mode)for(bool twice:{false,true})
    for(double overlap:{0.0,0.2,-0.1})for(bool anchored:{false,true}){
      std::vector<FieldSegment> segments;
      if(mode!=2)segments=base;
      if(mode!=0){segments.push_back(rib);segments.push_back(rib2);}
      if(mode==3){segments[0].width=1e-6;segments.back().trim=0.8;}
      std::vector<zScalarArray> a,b;zScalarArray sa,sb;
      reference::computeSegmentDistances(grid,segments,a,&sa);
      current::computeSegmentDistances(grid,segments,b,&sb);
      assert(a==b && sa==sb);
      zScalarArray wa,wb,ia,ib,pa,pb,ta,tb;
      reference::computeWidthField(grid,a,segments,twice,overlap,&wa,&ia);
      current::computeWidthField(grid,b,segments,twice,overlap,&wb,&ib);
      equal(wa,wb);equal(ia,ib);
      zObjGraph ga,gb;
      reference::computePathField(grid,boundary,segments,a,twice,overlap,&ga,&pa,&ta,
                                  anchored?&slot:nullptr,0.8,&sa);
      current::computePathField(grid,boundary,segments,b,twice,overlap,&gb,&pb,&tb,
                                anchored?&slot:nullptr,0.8,&sb);
      equal(pa,pb);equal(ta,tb);equalGraph(ga,gb);
      zScalarArray only;
      current::computeWidthField(grid,b,segments,twice,overlap,&only,nullptr);equal(only,wb);
      current::computeWidthField(grid,b,segments,twice,overlap,nullptr,&only);equal(only,ib);
      current::computePathField(grid,boundary,segments,b,twice,overlap,nullptr,&only,nullptr,
                                anchored?&slot:nullptr,0.8,nullptr);equal(only,pb);
      current::computePathField(grid,boundary,segments,b,twice,overlap,nullptr,nullptr,&only,
                                anchored?&slot:nullptr,0.8,&sb);equal(only,tb);
      comparisons+=int(grid.points.size());
    }
  zScalarArray path,track;zObjGraph slots=slot;
  current::computePathField(grid,boundary,{}, {},true,0,&slots,&path,&track,nullptr,0,nullptr);
  assert(slots.points.empty() && slots.edges.empty());
  for(double v:path)assert(v==1);for(double v:track)assert(v==0);
  std::vector<zScalarArray> shortDistances(4,zScalarArray(1));
  current::computePathField(grid,boundary,base,shortDistances,true,0,nullptr,&path,&track,nullptr,0,nullptr);
  for(double v:path)assert(v==1);
  // Reused diagnostics must clear when a later call produces no slots.
  slots=slot;std::vector<zScalarArray> distances;
  current::computeSegmentDistances(grid,base,distances,nullptr);
  current::computePathField(grid,boundary,base,distances,false,0,&slots,&path,&track,nullptr,0,nullptr);
  assert(slots.points.empty());
  for(bool ccw:{false,true})for(double offset:{-20.0,-1.0,0.0,1.0,20.0}){
    LoopAnchor anchor;anchor.edge=1;anchor.t=0.25;anchor.ccw=ccw;
    zObjGraph a,b;zPoint ca,cb;
    assert(reference::createAnchoredSlot(boundary.points,anchor,offset,1,a,&ca));
    assert(current::createAnchoredSlot(boundary.points,anchor,offset,1,b,&cb));
    equalGraph(a,b);assert(vecDistSq(ca,cb)<1e-20);
  }
  zPointArray tiny={{0,0},{1e-8,0},{1e-8,1e-8},{0,1e-8}};
  LoopAnchor anchor;anchor.edge=0;anchor.t=0.3;
  zObjGraph a,b;
  assert(reference::createAnchoredSlot(tiny,anchor,1e-9,1,a,nullptr));
  assert(current::createAnchoredSlot(tiny,anchor,1e-9,1,b,nullptr));equalGraph(a,b);
  std::cout<<"PASS: "<<comparisons<<" field-point comparisons, optional outputs, empty/short inputs, "
           <<"stale slots and wrapped/tiny anchored slots\n";
}
'''
with tempfile.TemporaryDirectory(prefix="ts3dp_composition_") as tmp:
    cpp = Path(tmp) / "field_composition.cpp"
    exe = Path(tmp) / "field_composition.exe"
    cpp.write_text(prelude + stubs + shared + "namespace reference {\n" + implementations(old)
                   + "}\nnamespace current {\n" + implementations(new) + "}\n" + checks)
    subprocess.run(["cl", "/nologo", "/EHsc", "/std:c++17", "/fp:precise", str(cpp),
                    "/Fe:" + str(exe), "/Fo:" + str(Path(tmp)/"field_composition.obj")],
                   check=True, cwd=tmp)
    subprocess.run([str(exe)], check=True)
