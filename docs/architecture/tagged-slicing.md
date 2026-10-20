# Tagged slicing and PrintSyn

The tagged slicing engine is available through `<zspace/slicing.h>` as
`zSpace::zTsSlicing`. It ports the topology/harmonic slicers and PrintSyn from
the desktop `slicer` branch (`2db24962`) into the maintained v3 source layout.
The original v3 `zTs3DP`, `z3DPSlicer` and `z3DPPrintSynthesis` APIs are unchanged.
The new name avoids colliding with that existing, differently shaped `zTs3DP` API.

Public headers live in `include/zspace/zToolsets/fabrication/slicing/` and
implementations in `src/zToolsets/fabrication/slicing/`. Headers no longer include
implementation `.cpp` files; link against `zspace::toolsets`.

## Input contract

- A polygon mesh, copied into the session. Topology mode requires all-quad faces,
  one connected floor patch and regular columns connecting the floor to the top.
- One `VertexTag` per **face**, in polygon order: `Bottom`, `Top`, or `Side`.
  Operators derive their vertex tags from these face tags.
- Optional per-face width-group IDs, indexing `SlicerSettings::BoundaryWidth`.
- An optional bracing graph in the input mesh's coordinate system.
- An explicit slot anchor in that coordinate system, required by both operators.
- Slicing and PrintSyn settings. Coordinates, layer spacing and widths use the
  same model units; the compute engine does not convert units.

`zTsSlicing` owns copies of the mesh, face tags, width IDs, bracing graph and
anchor. The lower-level operators borrow their inputs. Keep those inputs alive
when using `TopologySlicer` or `HarmonicSlicer` directly.

## Minimal usage

```cpp
#include <zspace/slicing.h>
#include <stdexcept>

zSpace::zTsSlicing session;
session.setInputMesh(mesh, faceTags, faceWidthIds);
session.setBracingGraph(bracing); // optional
session.setSlotAnchor(anchor);
session.getSettings().LayerSpacing = 1.0f;
session.getSettings().BoundaryWidth = {0.5f};
session.getSettings().FieldCellSize = 0.05f;
session.getPrintSynSetting().ResampleLength = 0.1f;

if (!session.computeSlicer(false))
    throw std::runtime_error(session.getStatus());
if (!session.computePostProcess(false))
    throw std::runtime_error(session.getStatus());

const auto& layers = session.get3DPrintPaths();
const auto& combined = session.getResultantPrintPath();
```

`computeBasicStructure()` validates/builds the initial structure.
`computeSlicer()` runs through layer meshes, unrolled charts and fields.
`computePrintPath()` extracts/maps paths without PrintSyn.
`computePostProcess()` extracts paths and applies PrintSyn resampling, seam
alignment, projection and linking. Failures retain diagnostic partial results
where the desktop implementation did so; check the return value and status
before consuming the combined path.

Per-path data includes the graph, widths, heights, normals and track IDs.
Topology intermediate outputs are available through `getTopologySlicer()`;
harmonic outputs through `getHarmonicSlicer()`. Changing mode replaces the
operator and invalidates references to its results.

## Build configurations

The default native build includes topology slicing and PrintSyn without TetGen,
OpenUSD or OpenMP. Enable tests with `ZSPACE_TOOLSETS_BUILD_TESTS=ON`.

```powershell
cmake -S . -B build/slicing -DZSPACE_TOOLSETS_BUILD_TESTS=ON
cmake --build build/slicing --config Release --parallel 4
ctest --test-dir build/slicing -C Release --output-on-failure
```

Run from a compiler-enabled shell, using the same generator/toolchain as the
rest of your v3 build. The native smoke suite includes topology/PrintSyn,
input validation, width groups/bracing/boundary slots, serial/parallel parity,
repeat computation, and moved-session reuse, alongside the
existing v3 3DP regression.

For harmonic slicing, build/provide TetGen for the same platform and runtime:

```powershell
cmake -S . -B build/slicing -DZSPACE_TOOLSETS_WITH_TETGEN=ON `
  -DZSPACE_TETGEN_ROOT=C:/SDKs/tetgen
```

The root should contain `include/tetgen.h` and `lib/tetgen.lib` (or the platform
equivalent). Separate installations can set `ZSPACE_TETGEN_INCLUDE_DIR` and
`ZSPACE_TETGEN_LIBRARY` explicitly. Nothing is downloaded or vendored by this
option. Harmonic sources and their regression are included only when enabled.
Query `zTsSlicing::supportsSlicerMode()` before selecting a mode. Selecting an
unavailable or invalid mode throws `std::invalid_argument` and preserves the
current session.

`ZSPACE_TOOLSETS_SLICING_OPENMP=ON` enables existing parallel loops. Sequential
calls (`useParallel=false`) remain available in either build. Eigen stays sequential
to avoid nested parallel scheduling.

For the optional desktop USD adapter, set
`ZSPACE_TOOLSETS_SLICING_WITH_OPENUSD=ON` and `ZSPACE_OPENUSD_ROOT` to an external
OpenUSD SDK. The adapter reuses core's OpenUSD CMake helper. Deploy the SDK's
runtime libraries and schema/plugin resources with the application. Its public
API is `zTsIO::readSlicerInputUSD()` / `exportPrintPathUSD()`; without the option,
these functions return a diagnostic explaining the missing capability.

The USD adapter preserves `/block`, `/bracing`, `/slotpoint`, face/width tags,
settings and print-path attributes. It is separate from the compute engine.

## Browser boundary

The dedicated `wasm/slicing` target compiles the same engine and PrintSyn sources
into `wasm/out/zspace_slicing.js` and `.wasm`. It leaves the established core and
toolsets runtimes unchanged. Run `wasm/scripts/build-slicing-wasm.ps1` with the
Emscripten SDK installed. The default TetGen source directory is
`wasm/cache/tetgen`; provide `-TetgenDir` for another source checkout. The reviewed
build uses Emscripten 6.0.13 and libigl/tetgen commit
`e05aca7df74e3f531bc35733ed87d36d437266c5`. Without TetGen, capability bit 2 is
absent and harmonic mode is rejected. OpenMP is disabled in this worker build.

The exported ABI has version 1 and capability bits 1 (topology) / 2 (harmonic).
Create/destroy manage handles; `zspace_slicing_compute_from_buffers` uploads
owned mesh arrays, face tags, optional width groups and bracing, anchor, mode and
settings using the existing ZST1 binary codec. Read results with
`zspace_slicing_result_buffers` and `zspace_toolsets_result_buffer_size`. Input
allocation uses `zspace_alloc`/`zspace_free`. Errors use
`zspace_toolsets_last_error`. Failed computation preserves the previous result.

Results contain owned slice/chart meshes, path/width fields, per-layer print
paths and the combined PrintSyn path with widths/heights/normals/track IDs. The
viewer computes in a module worker and copies the result before releasing its
native session. No live heap views reach viewer scene objects.

The separate `wasm/slicing-usd` target adds local USD input through the native
`readSlicerInputUSD` reader. Build the core OpenUSD WASM SDK first, then run
`wasm/scripts/build-slicing-usd-wasm.ps1`. This emits its own pair under
`wasm/out/slicing-usd/`. Tagged input preserves `/block`, `/bracing`, `/slotpoint`,
face/width tags and both settings dictionaries. Plain USD accepts `/block` or a
single polygon mesh, bakes transforms and rotates Y-up geometry to Z-up. Stage
units are preserved. Use self-contained USD or USDZ for referenced assets. The
importer uses two pthread workers and requires browser isolation headers.

In `zspace_alice_webviewer`, install both matching generated pairs under `wasm/`
and `wasm/slicing-usd/`, then select **WASM Slicing**. **Load local OBJ / USD**
opens a browser file picker; plain inputs infer cap tags from horizontal Z
faces. Controls support topology/harmonic computation, bracing and display of
layers, fields and combined paths. The viewer setup and installed-runtime tests
are documented in `docs/development/slicing.md` there. Generated pairs remain
ignored; source is maintained on the `slicer` feature branch for review.
### Legacy USD input compatibility

The reader also accepts custom `FaceTag` / `widthid` arrays when the equivalent
uniform primvars are absent. Authored primvars take precedence and invalid
primvars are rejected. Legacy Mesh bracing stores `edgeVertexIndices` pairs
with one width ID per edge; each unbranched connected chain/loop must use one
width ID. A Mesh containing one point is accepted for `slotpoint` as well.
The reader keeps authored units and all slicer/PrintSyn settings. The viewer's
length controls use the native limits, preserving millimetre widths above 10.