"""Run the actual multi-mesh transfer helper with OpenMP enabled.

Geometry adapters record the OpenMP team size and model a deterministic
barycentric transfer, verifying execution mode and identical mapped output.
"""
from pathlib import Path
import os
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
vcvars = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)')) / 'Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat'
setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True, check=True, capture_output=True, text=True)
for line in setup.stdout.splitlines():
    key, sep, value = line.partition('=')
    if sep and key:
        os.environ[key] = value

source = (ROOT/'cpp/source/zToolsets/ts3dp/utils/Utils.cpp').read_text()
kernel = source[source.index('ZSPACE_INLINE bool mapPointToMeshes('):source.index('ZSPACE_INLINE bool mapPoint(')]
prelude = r'''
#include <omp.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#define ZSPACE_INLINE
struct zPoint{double x=0,y=0,z=0;};
using zPointArray=std::vector<zPoint>;
using zIntArray=std::vector<int>;
struct zObjGraph{zPointArray points;};
struct zObjMesh{int id=-1;double scale=1;};
struct ParameterGraph{zPoint point;};
struct zFnGraph{
  zObjGraph &graph;
  explicit zFnGraph(zObjGraph &g):graph(g){}
  void create(zPointArray &p,zIntArray &){graph.points=p;}
  void getVertexPositions(zPointArray &p){p=graph.points;}
};
std::vector<int> teams(256);
bool parameterizeGraph(const zObjGraph &g,const zObjMesh &,ParameterGraph &out){out.point=g.points[0];return true;}
bool realizeParameterGraph(const ParameterGraph &p,zObjMesh &m,zObjGraph &out){
  if(m.id>=0)teams[m.id]=omp_get_num_threads();
  out.points={{p.point.x*m.scale,p.point.y*m.scale,p.point.z*m.scale}};return true;
}
double vecDistSq(const zPoint &a,const zPoint &b){double x=a.x-b.x,y=a.y-b.y,z=a.z-b.z;return x*x+y*y+z*z;}
'''
checks = r'''
int main(){
  omp_set_dynamic(0);omp_set_num_threads(4);
  zObjMesh source;
  std::vector<zObjMesh> meshes(256);std::vector<const zObjMesh*> targets(256);
  for(int i=0;i<256;++i){meshes[i]={i,double(i+1)};targets[i]=&meshes[i];}
  targets[17]=nullptr;
  zPoint point={0.2,0.3,0.4};zPointArray serial,parallel;std::vector<char> sv,pv;
  double distance=-1;
  assert(mapPointToMeshes(point,source,targets,serial,sv,&distance,false));
  assert(distance==0);
  for(int i=0;i<256;++i)if(i!=17)assert(teams[i]==1);
  std::fill(teams.begin(),teams.end(),0);
  assert(mapPointToMeshes(point,source,targets,parallel,pv,&distance,true));
  assert(sv==pv&&sv[17]==0);
  for(int i=0;i<256;++i){assert(vecDistSq(serial[i],parallel[i])==0);if(i!=17)assert(teams[i]==4);}
  std::cout<<"PASS: false uses 1 thread, true uses 4 threads; 256 targets have identical output\n";
}
'''
with tempfile.TemporaryDirectory(prefix='parallel_point_mapping_') as tmp:
    tmp = Path(tmp)
    cpp = tmp/'test.cpp'
    cpp.write_text(prelude + kernel + checks)
    exe = tmp/'test.exe'
    subprocess.run(['cl','/nologo','/EHsc','/std:c++17','/openmp',str(cpp),'/Fe:'+str(exe),'/Fo:'+str(tmp/'test.obj')],check=True,cwd=tmp)
    env = dict(os.environ, OMP_THREAD_LIMIT='4')
    subprocess.run([str(exe)],check=True,env=env)
