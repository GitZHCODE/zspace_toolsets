"""Syntax-check changed slicer/PrintSyn translation units with the local SDK.

Uses MSVC Build Tools and the adjacent zspace_alice2 SDK checkout. This is a
real-header compile check, not a link or end-to-end geometry pipeline test.
"""
from pathlib import Path
import os, subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--openmp', action='store_true', help='Check OpenMP pragmas with MSVC /openmp')
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
vcvars = Path('C:/Program Files (x86)/Microsoft Visual Studio/2022/BuildTools/VC/Auxiliary/Build/vcvars64.bat')
setup = subprocess.run(f'call "{vcvars}" >nul && set', shell=True, check=True, capture_output=True, text=True)
for line in setup.stdout.splitlines():
    k, sep, v = line.partition('=')
    if sep and k: os.environ[k] = v
sdk = root / 'zspace_alice2/alice2/sdk/zspace/include'
includes = [root/'zspace_toolsets/cpp', sdk, sdk/'depends',
            root/'zspace_alice2/alice2/depends/eigen',
            root/'zspace_toolsets/cpp/depends/igl/headers',
            root/'zspace_toolsets/cpp/depends/igl/build/_deps/tetgen-src']
sources = ['zSlicer/harmonic/'+n+'.cpp' for n in ['Structure','Heat','Harmonic','TracingUV','Contour','Field','PostProcess']]
sources += ['zSlicer/topology/'+n+'.cpp' for n in ['Structure','Unroll','Field','PostProcess']]
sources += ['printsyn/zPrintSyn.cpp','utils/Utils.cpp','zTs3DP.cpp']
for src in sources:
    print('CHECK',src,flush=True)
    cmd = ['cl','/nologo','/Zs','/EHsc','/std:c++17','/DZSPACE_TOOLSETS_STATIC_LIBRARY','/DZSPACE_STATIC_LIBRARY']
    if args.openmp:
        cmd.append('/openmp')
    cmd += ['/I'+str(p) for p in includes]
    cmd += [str(root/'zspace_toolsets/cpp/source/zToolsets/ts3dp'/src)]
    result = subprocess.run(cmd,cwd=Path(__file__).resolve().parent,capture_output=True,text=True)
    if result.returncode:
        print(result.stdout[-16000:]); print(result.stderr); raise SystemExit(result.returncode)
print('PASS: real SDK syntax checks (OpenMP '+('enabled' if args.openmp else 'disabled')+')')
