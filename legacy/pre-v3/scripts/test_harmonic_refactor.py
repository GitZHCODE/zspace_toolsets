"""Compile actual old/new harmonic kernels and compare their results.

Usage: py -3 scripts/test_harmonic_refactor.py PATH_TO_BEFORE_DIRECTORY
The directory contains HarmonicSlicer.h, Harmonic.cpp, Heat.cpp and Contour.cpp from
before the refactor. Eigen is real; mesh creation is a lightweight test adapter.
"""
from pathlib import Path
import os
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
BEFORE = Path(sys.argv[1])
CURRENT = ROOT / 'cpp/source/zToolsets/ts3dp/zSlicer/harmonic'
HEADER = ROOT / 'cpp/headers/zToolsets/ts3dp/zSlicer/HarmonicSlicer.h'
vcvars = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)')) / 'Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat'
setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True, check=True, capture_output=True, text=True)
for line in setup.stdout.splitlines():
    k, sep, v = line.partition('=')
    if sep and k:
        os.environ[k] = v

prelude = r'''
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#define ZSPACE_INLINE
struct zPoint {
  double x=0,y=0,z=0;
  zPoint()=default;
  zPoint(double a,double b,double c=0):x(a),y(b),z(c){}
  zPoint operator+(const zPoint &b)const{return {x+b.x,y+b.y,z+b.z};}
  zPoint operator-(const zPoint &b)const{return {x-b.x,y-b.y,z-b.z};}
  zPoint operator*(double b)const{return {x*b,y*b,z*b};}
  double operator*(const zPoint &b)const{return x*b.x+y*b.y+z*b.z;}
  zPoint operator^(const zPoint &b)const{return {y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x};}
  zPoint &operator+=(const zPoint &b){x+=b.x;y+=b.y;z+=b.z;return *this;}
};
using zVector=zPoint;
using zPointArray=std::vector<zPoint>;
using zIntArray=std::vector<int>;
using zFloatArray=std::vector<float>;
struct zObjMesh {zPointArray positions;zIntArray counts,connects;};
using zObjMeshArray=std::vector<zObjMesh>;
struct zFnMesh {
  zObjMesh &mesh;
  explicit zFnMesh(zObjMesh &m):mesh(m){}
  void clear(){mesh={};}
  void create(zPointArray &p,zIntArray &n,zIntArray &c){mesh={p,n,c};}
};
struct TriangleKeyHash {size_t operator()(const std::array<int,3> &a)const{return a[0]*73856093u^a[1]*19349663u^a[2]*83492791u;}};
unsigned long long makeEdgeKey(int a,int b){if(a>b)std::swap(a,b);return (static_cast<unsigned long long>(a)<<32)|unsigned(b);}
double triangleArea2(const zPoint &a,const zPoint &b,const zPoint &c){auto n=(b-a)^(c-a);return std::sqrt(n*n);}
void printError(const std::string &){}
enum class VertexTag {Bottom=0,Top=1,Side=2};
struct Settings {int _HEqualiseIterations=5;};
struct TriangleSoup {zPointArray Positions;std::vector<std::array<int,3>> Triangles;};
class zTetrahedron {
public:zPointArray Vertices;std::vector<std::array<int,4>> Tets;zIntArray Markers,BottomVertexIds,TopVertexIds;};
struct TetElement {bool Valid=false;double Volume=0;Eigen::Matrix<double,3,4> Gradients;Eigen::Matrix4d Stiffness;};
struct SurfaceOperators {Eigen::SparseMatrix<double> Laplacian;Eigen::VectorXd Mass;double MeanEdgeLength=0;};
struct PrintPath{int value=0;};
using PrintPathArray=std::vector<PrintPath>;
struct zObjMeshScalarField{int size=0;};
struct zFnMeshScalarField{zObjMeshScalarField &f;explicit zFnMeshScalarField(zObjMeshScalarField &v):f(v){}void clear(){f.size=0;}};
Eigen::Vector3d toEigenVec(const zPoint &p){return {p.x,p.y,p.z};}
double meanEdgeLength(const zPointArray &points,const std::vector<std::array<int,3>> &triangles){
  double sum=0;int n=0;std::unordered_set<unsigned long long> seen;
  for(auto tri:triangles)for(int i=0;i<3;++i){int a=tri[i],b=tri[(i+1)%3];if(!seen.insert(makeEdgeKey(a,b)).second)continue;
    auto d=points[a]-points[b];sum+=std::sqrt(d*d);++n;}
  return n?sum/n:0;
}
'''

def kernels(directory, header, namespace, old):
    hd = header.read_text()
    solver = hd[hd.index('  bool solveTetDirichlet(') if old else hd.index('  struct TetDirichletSystem {'):hd.index('  bool buildBottomParameterization();')]
    solver += hd[hd.index('  static bool buildSurfaceOperators('):hd.index('  bool solveTetDirichlet(') if old else hd.index('  struct TetDirichletSystem {')]
    contour_end = '  void extractHarmonicIsoMeshes(' if '  void extractHarmonicIsoMeshes(' in hd else '  struct WallWidthLookup {'
    contour = hd[hd.index('  int crossingVertex(') if old else hd.index('  struct SliceData {'):hd.index(contour_end)]
    contour = contour.replace('zTsUtils::TriangleKeyHash', 'TriangleKeyHash')
    harmonic = (directory / 'Harmonic.cpp').read_text()
    harmonic = harmonic[harmonic.index('ZSPACE_INLINE bool HarmonicSlicer::'):harmonic.index('// The harmonic field w')]
    iso = (directory / 'Contour.cpp').read_text()
    iso_end = 'ZSPACE_INLINE void HarmonicSlicer::extractHarmonicIsoMeshes(' if 'ZSPACE_INLINE void HarmonicSlicer::extractHarmonicIsoMeshes(' in iso else 'namespace {\n\nvoid fillUnresolvedIds'
    iso = iso[iso.index('ZSPACE_INLINE int HarmonicSlicer::'):iso.index(iso_end)]
    heat = (directory / 'Heat.cpp').read_text()
    heat = heat[heat.index('ZSPACE_INLINE bool HarmonicSlicer::buildSurfaceOperators'):heat.index('// d_b + d_t')]
    reset_members = reset_declarations = reset = ''
    if not old:
        reset_members = '''
zPointArray _BottomUV; zPoint _BottomUVScale, _BottomUVOrigin;
zIntArray o_TraceGraphs, _SurfaceTriWidthIds, i_VertexTags, _HarmonicIsoEdgeWidthIds;
std::pair<zFloatArray,zFloatArray> _HeatValue;
zObjMeshArray o_IsoMeshes,o_IsoMeshes2D;
std::vector<zPointArray> o_IsoMeshUV;
zIntArray o_IsoMeshBoundary,o_IsoMeshBoundary2D,o_IsoMeshBoundaryWidthId,o_BracingGraphs2D;
std::vector<std::vector<double>> o_PathFields,o_WidthFields,o_TrackFields;
zObjMeshScalarField o_ScalarField;
PrintPathArray o_PrintPaths2D,o_PrintPaths3D;
PrintPath o_ResultantPrintPath;
'''
        reset_declarations = hd[hd.index('  void clearHeat();'):hd.index('  void buildSurfaceSoup();')]
        structure = (directory/'Structure.cpp').read_text()
        reset = structure[structure.index('ZSPACE_INLINE void HarmonicSlicer::clearHeat()'):structure.index('//--------------------------',structure.index('ZSPACE_INLINE void HarmonicSlicer::clearHeat()'))]
        reset += structure[structure.index('ZSPACE_INLINE void HarmonicSlicer::clearPrintPaths()'):structure.index('ZSPACE_INLINE bool HarmonicSlicer::computeBasicStructure(')]
    return f'namespace {namespace} {{\nclass HarmonicSlicer {{public:\nzTetrahedron _Tet; std::vector<TetElement> _TetElements; Settings i_Settings; TriangleSoup _SurfaceSoup; zPointArray _ParamUVW; std::vector<bool> _ParamValid;\n' + reset_members + reset_declarations + solver + contour + '};\n' + harmonic + iso + heat + reset + '\n}\n'

checks = r'''
template<class S> void fill(S &s,std::mt19937 &rng) {
  std::uniform_real_distribution<double> d(0.2,2.0);
  s._Tet.Vertices={{0,0,0},{d(rng),0,0},{0,d(rng),0},{0,0,d(rng)},{0.1,0.1,0.1}};
  s._Tet.Tets={{{0,1,2,4}},{{0,1,4,3}},{{0,4,2,3}},{{4,1,2,3}}};
  s._Tet.BottomVertexIds={0,1,2};s._Tet.TopVertexIds={3};
  s._Tet.Markers={0,0,0,1,2};
  for(auto tet:s._Tet.Tets) {
    Eigen::Matrix4d m;
    for(int j=0;j<4;++j){auto p=s._Tet.Vertices[tet[j]];m.row(j)<<1,p.x,p.y,p.z;}
    TetElement e;e.Valid=true;e.Gradients=m.inverse().bottomRows<3>();e.Volume=std::abs(m.determinant())/6.;
    e.Stiffness=e.Volume*e.Gradients.transpose()*e.Gradients;
    s._TetElements.push_back(e);
  }
}
void compareField(const zFloatArray &a,const zFloatArray &b) {
  assert(a.size()==b.size());for(size_t i=0;i<a.size();++i)assert(std::abs(a[i]-b[i])<2e-6);
}
void comparePoints(const zPointArray &a,const zPointArray &b){
  assert(a.size()==b.size());for(size_t i=0;i<a.size();++i){auto d=a[i]-b[i];assert(d*d<1e-20);}
}
int main(){
  std::mt19937 rng(231);
  std::uniform_real_distribution<double> d(0,1);
  for(int c=0;c<200;++c){
    before::HarmonicSlicer old;fill(old,rng);
    after::HarmonicSlicer now;now._Tet=old._Tet;now._TetElements=old._TetElements;
    TriangleSoup soup;soup.Positions=old._Tet.Vertices;soup.Positions.pop_back();
    soup.Triangles={{{0,1,2}},{{0,3,1}},{{0,2,3}},{{1,3,2}}};
    SurfaceOperators ops;int skipped=0;assert(old.buildSurfaceOperators(soup,ops,skipped));
    Eigen::SparseMatrix<double> heatMatrix=ops.Laplacian*(ops.MeanEdgeLength*ops.MeanEdgeLength);
    for(int i=0;i<4;++i)heatMatrix.coeffRef(i,i)+=ops.Mass[i];
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> heatSolver;heatSolver.compute(heatMatrix);
    std::string heatError;
    for(const auto &sources:std::vector<zIntArray>{{0,1,2},{3},{0,1,2,3}}){
      zFloatArray a,b;assert(old.computeHeatGeodesics(soup,ops,sources,a,heatError));
      assert(now.computeHeatGeodesics(soup,ops,heatSolver,sources,b,heatError));compareField(a,b);
    }
    auto fixed=old.capDirichletValues();
    if(c%3==0){fixed[1]=-1;fixed[2]=-1;}
    if(c%7==0) fixed[4]=0.25;
    after::HarmonicSlicer::TetDirichletSystem system;
    std::string error;assert(now.prepareTetDirichlet(fixed,system,error));
    zFloatArray a,b;int free=0;
    assert(old.solveTetDirichlet(fixed,a,free,error));assert(now.solveTetDirichlet(system,b,error));compareField(a,b);
    std::vector<Eigen::Vector3d> targets(4);
    for(int k=0;k<10;++k){
      for(auto &t:targets)t=Eigen::Vector3d(d(rng),d(rng),d(rng));
      assert(old.solveTetDirichlet(fixed,a,free,error,&targets));assert(now.solveTetDirichlet(system,b,error,&targets));compareField(a,b);
    }
    // Equalisation uses fixed cap values, so prepare that exact system.
    fixed=old.capDirichletValues();assert(now.prepareTetDirichlet(fixed,system,error));
    assert(old.solveTetDirichlet(fixed,a,free,error));assert(now.solveTetDirichlet(system,b,error));
    assert(old.equaliseField(a,error));assert(now.equaliseField(system,b,error));compareField(a,b);
    for(int i=0;i<5;++i){old._ParamUVW.push_back({d(rng),d(rng),d(rng)});}
    now._ParamUVW=old._ParamUVW;
    if(c%5==0){old._Tet.Tets.push_back(old._Tet.Tets[0]);now._Tet.Tets=old._Tet.Tets;}
    zFloatArray phi;for(auto p:old._ParamUVW)phi.push_back(float(p.z));
    for(float iso:{0.f,0.1f,0.5f,0.9f,1.f,float(phi[0])}){
      zObjMeshArray meshes;std::vector<zPointArray> uvs;std::vector<std::vector<std::pair<int,int>>> origins;
      old.extractIsoMeshes(phi,{iso},meshes,uvs,origins);
      zObjMesh mesh;zPointArray uv;std::vector<std::pair<int,int>> edges;
      now.extractIsoMesh(iso,mesh,uv,edges);
      assert(mesh.connects==meshes[0].connects&&mesh.counts==meshes[0].counts);
      comparePoints(mesh.positions,meshes[0].positions);comparePoints(uv,uvs[0]);assert(edges==origins[0]);
    }
  }
  // A rejected cap triangle must never leave unreferenced vertices.
  after::HarmonicSlicer cap;cap._Tet.Markers={0,0,0};cap._SurfaceSoup.Positions={{0,0,0},{1,0,0},{0,1,0}};
  cap._SurfaceSoup.Triangles={{{0,1,2}}};cap._ParamUVW=cap._SurfaceSoup.Positions;cap._ParamValid={true,true,false};
  zObjMesh mesh;zPointArray uv;std::vector<std::pair<int,int>> origins;
  cap.extractCapMesh(VertexTag::Bottom,mesh,uv,origins);assert(mesh.positions.empty()&&uv.empty()&&origins.empty());
  cap._ParamValid[2]=true;cap.extractCapMesh(VertexTag::Bottom,mesh,uv,origins);assert(mesh.positions.size()==3&&mesh.counts.size()==1);
  // Re-tracing after a successful pipeline must clear stale derived results
  // while preserving the w field; replacing the structure also drops w.
  cap.o_IsoMeshes={mesh};cap.o_IsoMeshUV={uv};cap.o_IsoMeshes2D={mesh};
  cap.o_IsoMeshBoundary={1};
  cap.o_PathFields={{1}};cap.o_WidthFields={{1}};cap.o_TrackFields={{1}};
  cap.o_ScalarField.size=10;cap.o_PrintPaths2D={{1}};cap.o_PrintPaths3D={{1}};cap.o_ResultantPrintPath.value=1;
  cap._ParamUVW={{0.2,0.3,0.5}};cap._ParamValid={true};cap.o_TraceGraphs={1};cap._BottomUV={{0,0,0}};
  cap.clearParameterization();
  assert(cap.o_IsoMeshes.empty()&&cap.o_IsoMeshUV.empty()&&cap.o_IsoMeshes2D.empty());
  assert(cap.o_IsoMeshBoundary.empty());
  assert(cap.o_PathFields.empty()&&cap.o_WidthFields.empty()&&cap.o_TrackFields.empty()&&cap.o_ScalarField.size==0);
  assert(cap.o_PrintPaths2D.empty()&&cap.o_PrintPaths3D.empty()&&cap.o_ResultantPrintPath.value==0);
  assert(cap._ParamUVW.size()==1&&cap._ParamUVW[0].z==0.5&&std::isnan(cap._ParamUVW[0].x));
  assert(cap._ParamValid.size()==1&&!cap._ParamValid[0]&&cap.o_TraceGraphs.empty()&&cap._BottomUV.empty());
  cap.clearStructure();assert(cap._ParamUVW.empty()&&cap._ParamValid.empty()&&cap._SurfaceSoup.Triangles.empty()&&cap._Tet.Markers.empty());
  std::cout<<"PASS: 2200 Dirichlet solves, 200 equalisation runs, 600 heat solves, 1200 iso meshes; cap validation and downstream invalidation\n";
}
'''

eigen = ROOT.parent / 'zspace_alice2/alice2/depends/eigen'
with tempfile.TemporaryDirectory(prefix='harmonic_refactor_') as tmp:
    tmp = Path(tmp)
    cpp = tmp / 'test.cpp'
    cpp.write_text(prelude + kernels(BEFORE, BEFORE/'HarmonicSlicer.h', 'before', True)
                   + kernels(CURRENT, HEADER, 'after', False) + checks)
    exe = tmp / 'test.exe'
    subprocess.run(['cl','/nologo','/EHsc','/std:c++17','/fp:precise','/I'+str(eigen),
                    str(cpp),'/Fe:'+str(exe),'/Fo:'+str(tmp/'test.obj')],check=True,cwd=tmp)
    subprocess.run([str(exe)],check=True)
