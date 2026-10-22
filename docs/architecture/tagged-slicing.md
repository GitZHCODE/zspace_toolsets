# Tagged slicing and PrintSyn

The tagged slicing engine is available through `<zspace/slicing.h>` as
`zSpace::zTsSlicing`. It ports the topology/harmonic slicers and PrintSyn from
the desktop `slicer` branch (`2db24962`) into the maintained v3 source layout.
The earlier v3 3DP prototype has been removed. `zTsSlicing` is the maintained C++ slicing API.

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

Harmonic heights use the local normal-spacing estimate `deltaW / |grad(w)|`,
where `deltaW = 1 / (numberOfLayers - 1)` and `w` is the final harmonic field
(including any field equalisation). Tetrahedron volumes weight the nodal average;
heights are then interpolated on tetrahedral edges and barycentrically onto each
layer's paths. This is a local linear estimate, not an exact distance to the
previous curved section or a bound enforced by `LayerSpacing`. Both cap layers
use a full adjacent iso interval, matching topology's nonzero cap-height contract.
Zero/nonfinite gradients or unsupported vertices fail contour generation instead
of silently producing zero heights. `HarmonicSlicer::getVertexLayerHeights()`
returns one table per iso mesh, in that mesh's vertex order. Raw paths and PrintSyn
projection use these tables; linked paths retain the heights, so the existing
height-domain and print-volume statistics work in harmonic mode. Settings/input
invalidation clears the tables and statistics. Target-domain optimization remains
a separate operation.

`zTsUtils::buildPrintMesh()` lofts one centered rectangle per print-path vertex,
using its width, height and normalized print normal. Each graph component and
layer stays separate; open paths receive end caps. Zero widths/heights retain
collapsed sections, and invalid attributes or connectivity clear the mesh and
return a diagnostic with the layer/path location.

USD export preserves `/TDP_gcode` and also writes `/TDP_printmesh` (Mesh, display
name `print mesh`) and `/TDP_printpath` (linear BasisCurves, display name `print
path`). Pass the processed per-layer paths through the optional `layerPaths`
argument of `zTsIO::exportPrintPathUSD()` to avoid lofting travel connectors
between layers. Without it, the supplied path alone is rendered. Geometry and
widths are converted from internal millimetres to the source stage units;
re-export replaces these output prims while preserving the input stage.

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

Browser USD IO is owned by the core OpenUSD runtime. Toolsets WASM excludes
`utils/IO.cpp` and accepts only owned geometry, tags, width groups and settings.
The viewer translates `/AG_block`, `/Slicer_bracing`, `/Slicer_slot`,
`/Slicer_data` and `/Slicer_setting` through generic core USD queries and converts
lengths to millimetres. See the viewer's `docs/development/slicing.md` for builds
and tests. Native desktop IO remains available independently.
