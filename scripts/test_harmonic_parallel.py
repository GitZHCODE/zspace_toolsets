"""Compare real harmonic kernels in serial and OpenMP modes on tet fixtures.

Uses actual Eigen and source-extracted kernels; mesh/graph creation uses
lightweight adapters. Optional argument: snapshot repository root from before
parallelization, for comparison against the previous serial implementation.
"""
from pathlib import Path
import ast
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BEFORE = Path(sys.argv[1]) if len(sys.argv) > 1 else None
tree = ast.parse((ROOT/'scripts/test_harmonic_refactor.py').read_text())
prelude = next(ast.literal_eval(node.value) for node in tree.body
               if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'prelude' for t in node.targets))
prelude = prelude.replace('struct Settings {int _HEqualiseIterations=5;};',
    'struct Settings {int _HEqualiseIterations=5; double _HTraceStepFraction=0.2; int _HTraceMaxSteps=10000; double LayerSpacing=0.125; int _HArapIterations=3; int _HUnrollStart=0;};')
prelude = prelude.replace('  void clear(){mesh={};}', '''  void clear(){mesh={};}
  int numVertices(){return int(mesh.positions.size());}
  int numPolygons(){return int(mesh.counts.size());}
  void getVertexPositions(zPointArray &p){p=mesh.positions;}
  void getPolygonData(zIntArray &connects,zIntArray &counts){connects=mesh.connects;counts=mesh.counts;}''')
prelude += r'''
#include <omp.h>
#include <map>
#include <set>
struct zObjGraph{zPointArray positions;zIntArray edges;};
using zObjGraphArray=std::vector<zObjGraph>;
struct zFnGraph{
  zObjGraph &graph;
  explicit zFnGraph(zObjGraph &g):graph(g){}
  int numVertices(){return int(graph.positions.size());}
  void create(zPointArray &p,zIntArray &e){graph.positions=p;graph.edges=e;}
};
void printSuccess(const std::string &){}
struct TriElement {bool Valid=false;double Area=0;Eigen::Matrix3d Gradients=Eigen::Matrix3d::Zero();};
using FaceAdjacency=std::unordered_map<std::array<int,3>,std::array<int,2>,TriangleKeyHash>;
using SurfaceEdgeAdjacency=std::unordered_map<unsigned long long,std::array<int,2>>;
using EdgeWidthIds=std::unordered_map<unsigned long long,int>;
namespace zTsUtils{int layersFromSpacing(double length,double spacing){return std::max(1,int(std::ceil(length/spacing)))+1;}}
'''

# Use the real ARAP solver and boundary loop extraction, alongside mesh adapters.
utils = ROOT/'cpp/source/zToolsets/ts3dp/utils'
boundary = (utils/'Utils.cpp').read_text()
boundary = boundary[boundary.index('ZSPACE_INLINE double signedArea2D'):boundary.index('//---- HalfEdge Mesh Topology Walk Utilities')]
unroll = (utils/'Unroll.cpp').read_text()
unroll = unroll[unroll.index('namespace {', unroll.index('//---- As-rigid-as-possible')):unroll.rindex('} // namespace zTsUnroll')]
prelude += '''
namespace zTsUtils {
Eigen::Vector3d toEigenVec(const zPoint &p){return ::toEigenVec(p);}
double vecDistSq(const zPoint &a,const zPoint &b){auto d=a-b;return d*d;}
''' + boundary + '''}
using zTsUtils::signedArea2D;
using zTsUtils::outerBoundaryLoop;
using zTsUtils::layersFromSpacing;
namespace zTsUnroll {
struct ArapStats {double stretchMax=0,stretchMean=0;int flipped=0;};
''' + unroll + '}\n'

def kernels(repo, namespace, old=False):
    base = repo/'cpp'
    h = (base/'headers/zToolsets/ts3dp/zSlicer/HarmonicSlicer.h').read_text()
    s = base/'source/zToolsets/ts3dp/zSlicer/harmonic'
    declarations = h[h.index('  bool buildBottomParameterization();'):h.index('\n};',h.index('  bool buildBottomParameterization();'))]
    declarations = declarations.replace('zTsUtils::TriangleKeyHash','TriangleKeyHash')
    trace_path = s/'TracingUV.cpp'
    if not trace_path.exists():
        trace_path = ROOT/'cpp/source/zToolsets/ts3dp/zSlicer/harmonic/TracingUV.cpp'
    trace = trace_path.read_text()
    trace_parallel = 'computeVolumetricParameterization(bool' in trace
    trace = trace[trace.index('ZSPACE_INLINE bool HarmonicSlicer::buildBottomParameterization()'):trace.rindex('} // namespace zSpace')]
    contour = (s/'Contour.cpp').read_text()
    has_unroll = 'HarmonicSlicer::computeUnroll(' in contour
    contour = contour[contour.index('ZSPACE_INLINE int HarmonicSlicer::crossingVertex('):contour.rindex('} // namespace zSpace')]
    structure = (s/'Structure.cpp').read_text()
    structure = structure[structure.index('ZSPACE_INLINE void HarmonicSlicer::buildTetElements('):structure.index('ZSPACE_INLINE void HarmonicSlicer::clearHeat()')]
    prepare_parallel = 'bool useParallel' in structure
    heat_path = s/'Heat.cpp'
    if not heat_path.exists():
        heat_path = ROOT/'cpp/source/zToolsets/ts3dp/zSlicer/harmonic/Heat.cpp'
    heat = heat_path.read_text()
    heat = heat[heat.index('ZSPACE_INLINE bool HarmonicSlicer::columnLengths('):heat.index('ZSPACE_INLINE bool HarmonicSlicer::computeHeatDistances()')]
    members = r'''
zTetrahedron _Tet;std::vector<TetElement> _TetElements;TriangleSoup _SurfaceSoup;
zIntArray _SurfaceTriWidthIds;std::pair<zFloatArray,zFloatArray> _HeatValue;
zPointArray _BottomUV,_ParamUVW;zPoint _BottomUVScale,_BottomUVOrigin;
std::vector<bool> _ParamValid;zObjGraphArray o_TraceGraphs;std::string _Status;Settings i_Settings;
zObjMeshArray o_IsoMeshes;std::vector<zPointArray> o_IsoMeshUV;std::vector<EdgeWidthIds> _HarmonicIsoEdgeWidthIds;
zObjMeshArray o_IsoMeshes2D;
zObjGraphArray o_IsoMeshBoundary,o_IsoMeshBoundary2D;
std::vector<zIntArray> o_IsoMeshBoundaryWidthId;
void clearUnroll(){o_IsoMeshes2D.clear();o_IsoMeshBoundary.clear();o_IsoMeshBoundary2D.clear();o_IsoMeshBoundaryWidthId.clear();}
void clearContours(){clearUnroll();o_IsoMeshes.clear();o_IsoMeshUV.clear();_HarmonicIsoEdgeWidthIds.clear();}
void clearParameterization(){
  o_TraceGraphs.clear();_BottomUV.clear();_ParamValid.assign(_ParamUVW.size(),false);
  for(auto &p:_ParamUVW){p.x=p.y=std::numeric_limits<float>::quiet_NaN();}
}
bool columnLengths(double &,double &)const;
'''
    flags = ', bool useParallel' if prepare_parallel else ''
    # Stable adapters expose the same test API across the before/after snapshots.
    methods = f'''
void buildTetElements(int &count{flags});
bool computeVolumetricParameterization({'bool useParallel' if trace_parallel else ''});
bool computeContours(bool useParallel);
{'bool computeUnroll(bool useParallel);' if has_unroll else ''}
void prepare(bool parallel,int &count){{buildTetElements(count{',parallel' if prepare_parallel else ''});}}
bool trace(bool parallel){{return computeVolumetricParameterization({'parallel' if trace_parallel else ''});}}
void extract(bool parallel){{assert(computeContours(parallel));}}
'''
    return f'namespace {namespace} {{class HarmonicSlicer {{public:\n'+members+declarations+methods+'};\n'+structure+heat+trace+contour+'}\n'

checks = r'''
template<class S> void fixture(S &s,bool warped,int maxSteps){
  constexpr int n=4;
  auto index=[](int i,int j,int k){return (k*5+j)*5+i;};
  zPointArray points;
  for(int k=0;k<=n;++k)for(int j=0;j<=n;++j)for(int i=0;i<=n;++i){
    double x=i*0.5,y=j*0.25,z=k*0.25;
    if(warped)z*=1+0.15*x+0.1*y;
    points.emplace_back(x,y,z);
  }
  std::vector<std::array<int,4>> tets;
  for(int k=0;k<n;++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
    std::array<int,3> order={0,1,2};
    do{int c[3]={i,j,k};std::array<int,4> tet;tet[0]=index(c[0],c[1],c[2]);
      for(int step=0;step<3;++step){++c[order[step]];tet[step+1]=index(c[0],c[1],c[2]);}
      tets.push_back(tet);
    }while(std::next_permutation(order.begin(),order.end()));
  }
  std::map<std::array<int,3>,int> faceCounts;
  for(auto tet:tets)for(int omit=0;omit<4;++omit){std::array<int,3> face;int p=0;
    for(int q=0;q<4;++q)if(q!=omit)face[p++]=tet[q];std::sort(face.begin(),face.end());++faceCounts[face];}
  std::set<int> boundary;
  for(auto entry:faceCounts)if(entry.second==1)for(int id:entry.first)boundary.insert(id);
  std::vector<int> ordered(boundary.begin(),boundary.end()),newId(points.size());
  size_t surfaceCount=ordered.size();
  for(int id=0;id<int(points.size());++id)if(!boundary.count(id))ordered.push_back(id);
  for(int id=0;id<int(ordered.size());++id){
    int original=ordered[id];newId[original]=id;s._Tet.Vertices.push_back(points[original]);
    int k=original/25;int marker=k==0?0:(k==n?1:(boundary.count(original)?2:-1));
    s._Tet.Markers.push_back(marker);s._ParamUVW.emplace_back(0,0,k*0.25);
    if(marker==0)s._Tet.BottomVertexIds.push_back(id);if(marker==1)s._Tet.TopVertexIds.push_back(id);
    if(id<int(surfaceCount)){s._SurfaceSoup.Positions.push_back(points[original]);s._HeatValue.first.push_back(k*0.25f);s._HeatValue.second.push_back(1-k*0.25f);}
  }
  for(auto tet:tets){for(int &id:tet)id=newId[id];s._Tet.Tets.push_back(tet);}
  // Include a degenerate element, disconnected from every non-degenerate face.
  s._Tet.Tets.push_back({0,0,0,0});
  for(auto entry:faceCounts)if(entry.second==1){auto tri=entry.first;for(int &id:tri)id=newId[id];
    s._SurfaceSoup.Triangles.push_back(tri);
    bool bottom=true,top=true;for(int id:tri){bottom&=s._Tet.Markers[id]==0;top&=s._Tet.Markers[id]==1;}
    s._SurfaceTriWidthIds.push_back(bottom||top?-1:0);
  }
  s.i_Settings._HTraceMaxSteps=maxSteps;
}
bool sameNumber(double a,double b){return (std::isnan(a)&&std::isnan(b))||a==b;}
void pointsEqual(const zPointArray &a,const zPointArray &b){assert(a.size()==b.size());
  for(size_t i=0;i<a.size();++i)assert(sameNumber(a[i].x,b[i].x)&&sameNumber(a[i].y,b[i].y)&&sameNumber(a[i].z,b[i].z));}
template<class A,class B> void compare(const A &a,const B &b){
  assert(a._TetElements.size()==b._TetElements.size());
  for(size_t i=0;i<a._TetElements.size();++i){const auto &x=a._TetElements[i],&y=b._TetElements[i];
    assert(x.Valid==y.Valid&&x.Volume==y.Volume);if(x.Valid){assert(x.Gradients==y.Gradients);assert(x.Stiffness==y.Stiffness);}}
  assert(a._ParamValid==b._ParamValid);pointsEqual(a._ParamUVW,b._ParamUVW);
  assert(a.o_TraceGraphs.size()==b.o_TraceGraphs.size());
  for(size_t i=0;i<a.o_TraceGraphs.size();++i){pointsEqual(a.o_TraceGraphs[i].positions,b.o_TraceGraphs[i].positions);assert(a.o_TraceGraphs[i].edges==b.o_TraceGraphs[i].edges);}
  assert(a.o_IsoMeshes.size()==b.o_IsoMeshes.size());
  for(size_t k=0;k<a.o_IsoMeshes.size();++k){pointsEqual(a.o_IsoMeshes[k].positions,b.o_IsoMeshes[k].positions);
    assert(a.o_IsoMeshes[k].connects==b.o_IsoMeshes[k].connects&&a.o_IsoMeshes[k].counts==b.o_IsoMeshes[k].counts);
    pointsEqual(a.o_IsoMeshUV[k],b.o_IsoMeshUV[k]);assert(a._HarmonicIsoEdgeWidthIds[k]==b._HarmonicIsoEdgeWidthIds[k]);
    pointsEqual(a.o_IsoMeshes2D[k].positions,b.o_IsoMeshes2D[k].positions);
    assert(a.o_IsoMeshes2D[k].connects==b.o_IsoMeshes2D[k].connects&&a.o_IsoMeshes2D[k].counts==b.o_IsoMeshes2D[k].counts);
    pointsEqual(a.o_IsoMeshBoundary2D[k].positions,b.o_IsoMeshBoundary2D[k].positions);
    assert(a.o_IsoMeshBoundary2D[k].edges==b.o_IsoMeshBoundary2D[k].edges);
    assert(a.o_IsoMeshBoundaryWidthId[k]==b.o_IsoMeshBoundaryWidthId[k]);}
  for(size_t k=0;k<a.o_IsoMeshBoundary.size();++k){
    const auto &g=a.o_IsoMeshBoundary[k];assert(g.edges==a.o_IsoMeshBoundary2D[k].edges);
    assert(g.positions.size()==a.o_IsoMeshBoundary2D[k].positions.size());
    for(const auto &p:g.positions){bool found=false;for(const auto &q:a.o_IsoMeshes[k].positions)found|=p.x==q.x&&p.y==q.y&&p.z==q.z;assert(found);}
  }
}
int main(){
  omp_set_dynamic(0);
  for(bool warped:{false,true})for(int steps:{1,10000})for(int start:{0,1}){
    current::HarmonicSlicer serial;fixture(serial,warped,steps);int degenerate=0;
    serial.i_Settings._HUnrollStart=start;
    serial.prepare(false,degenerate);assert(degenerate==1);assert(serial.trace(false));
    std::string traceStatus=serial._Status;serial.extract(false);
    assert(serial.o_IsoMeshes.size()==9);
    size_t failed=std::count(serial._ParamValid.begin(),serial._ParamValid.end(),false);
    if(steps==1)assert(failed>0&&traceStatus.find("StepLimit")!=std::string::npos);
    if(!warped&&steps>1)assert(std::count(serial._ParamValid.begin(),serial._ParamValid.end(),true)>serial._Tet.BottomVertexIds.size());
    for(int threads:{1,2,4})for(int repeat=0;repeat<3;++repeat){
      omp_set_num_threads(threads);current::HarmonicSlicer parallel;fixture(parallel,warped,steps);
      parallel.i_Settings._HUnrollStart=start;
      parallel.prepare(true,degenerate);assert(degenerate==1);assert(parallel.trace(true));parallel.extract(true);compare(serial,parallel);
    }
    BEFORE_CHECK
  }
  std::cout<<"PASS: 72 serial/parallel comparisons at 1/2/4 threads; warped/flat fixtures, planar/UV ARAP, failed paths, tet matrices, 3D/2D iso meshes, boundary graphs and width IDs\n";
}
'''
vcvars = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)'))/'Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat'
setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True, check=True, capture_output=True, text=True)
for line in setup.stdout.splitlines():
    key, sep, value = line.partition('=')
    if sep and key:
        os.environ[key] = value
code = prelude + kernels(ROOT,'current')
if BEFORE:
    code += kernels(BEFORE,'before',True)
    checks = checks.replace('BEFORE_CHECK', 'before::HarmonicSlicer original;fixture(original,warped,steps);original.i_Settings._HUnrollStart=start;original.prepare(false,degenerate);assert(original.trace(false));original.extract(false);compare(serial,original);')
else:
    checks = checks.replace('BEFORE_CHECK','')
with tempfile.TemporaryDirectory(prefix='harmonic_parallel_') as tmp:
    tmp = Path(tmp)
    cpp = tmp/'test.cpp'
    cpp.write_text(code+checks)
    exe = tmp/'test.exe'
    eigen = ROOT.parent/'zspace_alice2/alice2/depends/eigen'
    subprocess.run(['cl','/nologo','/EHsc','/std:c++17','/openmp','/fp:precise','/I'+str(eigen),str(cpp),'/Fe:'+str(exe),'/Fo:'+str(tmp/'test.obj')],check=True,cwd=tmp)
    subprocess.run([str(exe)],check=True,env=dict(os.environ,OMP_THREAD_LIMIT='4'))
