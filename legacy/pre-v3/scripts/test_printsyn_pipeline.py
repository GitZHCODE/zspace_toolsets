"""Exercise the actual PrintSyn orchestrator with instrumented stage adapters."""
from pathlib import Path
import os
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
vcvars = Path(os.environ.get('ProgramFiles(x86)', 'C:/Program Files (x86)')) / 'Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat'
setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True, check=True, capture_output=True, text=True)
for line in setup.stdout.splitlines():
    k, sep, v = line.partition('=')
    if sep and k:
        os.environ[k] = v
source = (ROOT/'cpp/source/zToolsets/ts3dp/printsyn/zPrintSyn.cpp').read_text()
body = source[source.index('ZSPACE_INLINE bool zPrintSyn::computePipeline('):source.rindex('} // namespace zSpace')]
header = (ROOT/'cpp/headers/zToolsets/ts3dp/printsyn/zPrintSyn.h').read_text()
declarations = header[header.index('  enum class SeamMapping'):header.index('  /*! \\brief Step 1:')]
prelude = r'''
#include <cassert>
#include <iostream>
#include <string>
#include <vector>
#define ZSPACE_INLINE
using zScalarArray=std::vector<double>;
using zFloatArray=std::vector<float>;
struct PrintPath{int value=0;};
using PrintPathArray=std::vector<PrintPath>;
using zObjMeshArray=std::vector<int>;
struct zObjMeshScalarField{};
struct PrintSynSetting{double SeamOffset=0,SeamSize=0;};
int failStage=0,badStage=0,alignMode=0;
const zFloatArray *observedHeights=nullptr;
std::vector<int> calls;
std::vector<bool> flags;
bool stage(int n,PrintPathArray &out,std::string &status,std::vector<int>*bad,bool useParallel){
  flags.push_back(useParallel);calls.push_back(n);out={{n}};status="stage failure";
  if(bad){bad->clear();if(n==badStage)bad->push_back(2);}
  return n!=failStage;
}
class zPrintSyn{public:
'''
adapters = r'''
static bool extractPrintPaths(const zObjMeshScalarField &,const std::vector<const zScalarArray*>&,
  const std::vector<const zScalarArray*>&,const std::vector<const zScalarArray*>&,
  PrintPathArray &out,std::string &s,std::vector<int>*bad,bool useParallel=true){return stage(1,out,s,bad,useParallel);}
static bool resamplePrintPaths(PrintPathArray &out,const PrintSynSetting &,std::string &s,std::vector<int>*bad,bool useParallel=true){return stage(2,out,s,bad,useParallel);}
static bool alignStartPoints(PrintPathArray &out,const zObjMeshArray &,double,double,std::string &s,std::vector<int>*bad,bool useParallel=true){alignMode=1;return stage(3,out,s,bad,useParallel);}
static bool alignStartPoints(PrintPathArray &out,const zObjMeshArray &,const zObjMeshArray &,double,double,std::string &s,std::vector<int>*bad,bool useParallel=true){alignMode=2;return stage(3,out,s,bad,useParallel);}
static bool projectPrintPaths(const PrintPathArray &,const zObjMeshArray &,const zObjMeshArray &,
 const zFloatArray *h,PrintPathArray &out,std::string &s,std::vector<int>*bad=nullptr,bool useParallel=true){observedHeights=h;return stage(4,out,s,bad,useParallel);}
static bool linkPrintPaths(const PrintPathArray &,PrintPath &out,std::string &s){calls.push_back(5);out.value=5;s="link failure";return failStage!=5;}
};
'''
checks = r'''
int main(){
  zObjMeshScalarField grid;std::vector<zScalarArray> fields(3);zObjMeshArray charts(3);
  zFloatArray heights(3);PrintPathArray p2,p3;PrintPath combined;std::string status;
  PrintSynSetting settings;
  for(bool parallel:{false,true})for(auto mode:{zPrintSyn::SeamMapping::SharedTopology,zPrintSyn::SeamMapping::WorldNearest}){
    const zFloatArray *h=mode==zPrintSyn::SeamMapping::SharedTopology?&heights:nullptr;
    for(int failure=0;failure<=5;++failure)for(int bad=0;bad<=4;++bad){
      failStage=failure;badStage=bad;alignMode=0;calls.clear();flags.clear();observedHeights=reinterpret_cast<const zFloatArray*>(1);
      p2={{99}};p3={{99}};combined.value=99;
      bool ok=zPrintSyn::computePipeline(grid,fields,fields,fields,charts,charts,h,mode,settings,p2,p3,combined,status,parallel);
      int stop=failure?failure:6;if(bad&&bad<stop)stop=bad;
      assert(ok==(stop==6));
      std::vector<int> expected;for(int n=1;n<=std::min(stop,5);++n)expected.push_back(n);
      if(stop<4)expected.push_back(4);assert(calls==expected);for(bool flag:flags)assert(flag==parallel);
      assert(observedHeights==h);
      if(stop>=3)assert(alignMode==(mode==zPrintSyn::SeamMapping::SharedTopology?1:2));
      if(!ok){assert(combined.value==0);assert(!p2.empty());assert(status.find("stopped at ")==0);}
      else assert(combined.value==5);
      if(stop==4&&failure==4) assert(!p3.empty()); // preserve diagnostic projection output
      if(stop<4&&failure==4)assert(p3.empty()); // failed diagnostic projection clears stale 3D paths
    }
  }
  calls.clear();combined.value=99;p2={{99}};p3={{99}};
  assert(!zPrintSyn::computePipeline(grid,{},fields,fields,charts,charts,nullptr,zPrintSyn::SeamMapping::WorldNearest,settings,p2,p3,combined,status));
  assert(calls.empty()&&combined.value==0&&p2.empty()&&p3.empty());
  std::cout<<"PASS: both seam modes; stage failures, bad layers, partial results and invalid input (121 cases; serial/parallel flags including diagnostic projection)\n";
}
'''
with tempfile.TemporaryDirectory(prefix='printsyn_pipeline_') as tmp:
    tmp = Path(tmp)
    cpp = tmp/'test.cpp'
    cpp.write_text(prelude + declarations + adapters + body + checks)
    exe = tmp/'test.exe'
    subprocess.run(['cl','/nologo','/EHsc','/std:c++17',str(cpp),'/Fe:'+str(exe),'/Fo:'+str(tmp/'test.obj')],check=True,cwd=tmp)
    subprocess.run([str(exe)],check=True)
