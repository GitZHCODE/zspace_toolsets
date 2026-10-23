# Build the slicer engine

Use CMake 3.21 or newer and a C++17 toolchain. The default source build expects
`zspace_core` beside `zspace_toolsets`; set `ZSPACE_CORE_DIR` for another location.
These commands are documented build entry points, not results from a build run
performed for this documentation change.

## Native topology and PrintSyn

From a compiler-enabled shell in the toolsets repository:

```powershell
cmake -S . -B build/slicing -DZSPACE_TOOLSETS_BUILD_TESTS=ON
cmake --build build/slicing --config Release --parallel 4
ctest --test-dir build/slicing -C Release --output-on-failure
```

The default includes topology slicing and PrintSyn. Optional features are:

| CMake option | Requirement |
| --- | --- |
| `ZSPACE_TOOLSETS_WITH_TETGEN=ON` | Matching native TetGen headers/library; set `ZSPACE_TETGEN_ROOT`, or explicit include/library paths |
| `ZSPACE_TOOLSETS_SLICING_WITH_OPENUSD=ON` | External native OpenUSD SDK through `ZSPACE_OPENUSD_ROOT` |
| `ZSPACE_TOOLSETS_SLICING_OPENMP=ON` | Compiler OpenMP support |
| `ZSPACE_TOOLSETS_SLICING_ONLY=ON` | Restricts implementation sources to slicing |

Build libraries and SDKs for the same platform/runtime as the consumer. The
OpenUSD adapter also needs its runtime libraries and plugin/schema resources
when deployed. The TetGen option does not fetch an SDK automatically.

## Browser compute module

With Emscripten installed and a TetGen source checkout for harmonic support:

```powershell
& ./wasm/scripts/build-slicing-wasm.ps1 -CoreDir ../zspace_core -TetgenDir C:/SDKs/tetgen
```

The output pair is `wasm/out/zspace_slicing.js` and `zspace_slicing.wasm`.
Deploy both together. Without the expected TetGen sources, the build omits
harmonic support; capability bits identify topology (1) and harmonic (2).

At this source snapshot the compute target compiles with pthreads and OpenMP.
It needs SharedArrayBuffer in a secure, cross-origin-isolated browser context.
Serve COOP `same-origin` and COEP `require-corp` headers. The runtime sizes its
thread pool from browser-reported logical processors; this does not measure
current CPU load or guarantee a speedup. The compute binary excludes USD IO.

The web viewer snapshot has a **Compute slices** action, while the toolsets WASM
notes describe a newer **Parallel compute** action. That action is not present
in the reviewed web UI. Also, the web USD-importer build entry point referenced
by that viewer is absent from this toolsets commit. See the
[browser deployment guide](https://github.com/GitZHACODE/zspace_alice_webviewer/blob/feature/zspace-v3-gitbook/docs/gitbook/slicer/deployment.md)
before refreshing browser artifacts.

## Sources

- [Native CMake options and tests](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/CMakeLists.txt)
- [Slicing WASM target](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/wasm/slicing/CMakeLists.txt)
- [WASM build script](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/wasm/scripts/build-slicing-wasm.ps1)
