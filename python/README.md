# zSpace Python bindings

Core bindings live in `zspace_core/python/core`; slicer bindings live in `zspace_toolsets/python/toolsets/zTs3DP`. C++ algorithms remain in the v3
`include/zspace` and `src` trees. Python imports remain `zspace.core` and
`zspace.toolset.slicer`.

## Layout

- `build_headless.ps1`: Windows build entry point for Python 3.10 and 3.13.
- `CMakeLists.txt`: headless core, v3 toolsets and Python package configuration.
- `../../zspace_core/python/core/src`: geometry bindings maintained by Core.
- `toolsets/zTs3DP/src`: slicer bindings (historical directory name).
- `scripts/smoke_headless.py`: optional smoke example, not run automatically.
- `build_headless310/package` and `build_headless313/package`: generated packages.

## Build

Requires Visual Studio 2022 C++ tools, CMake 3.21+, x64 Python and pybind11.
The default core checkout is a sibling of `zspace_toolsets`:

```text
repo/
  zspace_core/
  zspace_toolsets/
    python/
```

From this `python` directory:

```powershell
python -m pip install pybind11 numpy
./build_headless.ps1 -PythonExecutable C:/Path/To/Python313/python.exe
```

Use Python 3.10 for the separate 3.10 package. Existing build caches from the old
repository must not be reused. `-BuildDirectory` and `-Pybind11Dir` override the
build location and pybind11 installation.

Topology is enabled by default. Add `-WithHarmonic -TetGenRoot C:/SDK/tetgen`
for Harmonic; that installation must provide `include/tetgen.h` and a compatible
library under `lib`. The script enables OpenMP. USD is disabled by default;
`-WithUSD` enables Core IO and the shared Toolsets slicing IO bindings.
The build copies SDK runtime DLLs and plugin resources beside `_core.pyd`.
Toolsets IO calls Core USDapi; SDK selection and deployment remain Core-owned.

CMake overrides include `ZSPACE_CORE_SOURCE`, `ZSPACE_TOOLSETS_SOURCE`,
`ZSPACE_PYTHON_BUILD_SLICER`, `ZSPACE_PYTHON_WITH_USD`,
`ZSPACE_TOOLSETS_WITH_TETGEN`, and `ZSPACE_TETGEN_ROOT`.
Core display, interop, tests and OpenGL are disabled. Core IO supplies the generic `USDdoc` binding. Without `-WithUSD`,
`USDdoc.supported()` returns false and USD operations report an error.

## Package and API compatibility

Set `PYTHONPATH` to the generated `package` directory. `_core.pyd` links native
core/interface/IO DLLs; `_slicer.pyd` links `zSpace_Toolsets.dll`.
These native libraries are copied into the package. Keep the package directories
together and use a Python interpreter matching the extension ABI.

The slicer class is `zTsSlicing`; `zTs3DP` remains a Python alias for existing
scripts. `zTsSlicing.supportsSlicerMode(mode)` reports whether a mode was compiled.
The old implementation-inclusion `impl.cpp` is removed: bindings now link the
same `zspace::toolsets` CMake target as native v3 clients.

## Bound API

The inventory below follows the current binding source, not every API available in the C++ SDK. Names retain their exact Python spelling. Overloads share one listed method name; see the linked binding source for argument names, constructors and return-value conversions. STL containers are converted through pybind11; modify and reassign list-valued properties when necessary rather than assuming a returned Python list is a live C++ container.

### zspace.core

#### `zVector`

Source: [bind_core_types.cpp](core/src/bind_core_types.cpp).

Properties: `x`, `y`, `z`.

Methods: `__eq__`, `__ne__`, `__getitem__`, `__add__`, `__sub__`, `__mul__`, `__xor__`, `__truediv__`, `__iadd__`, `__isub__`, `__imul__`, `__itruediv__`, `length2`, `length`, `normalize`, `squareDistanceTo`, `distanceTo`, `angle`, `angle360`, `dihedralAngle`, `cotan`, `getComponents`, `rotateAboutAxis`, `__repr__`.

#### `zColor`

Source: [bind_core_types.cpp](core/src/bind_core_types.cpp).

Properties: `r`, `g`, `b`, `a`, `h`, `s`, `v`.

Methods: `__eq__`, `toHSV`, `toRGB`, `__repr__`.

#### `zCurvature`

Source: [bind_core_types.cpp](core/src/bind_core_types.cpp).

Properties: `k1`, `k2`.

#### `zHEData`

Source: [bind_core_types.cpp](core/src/bind_core_types.cpp).

Enum values: `zVertexData`, `zHalfEdgeData`, `zEdgeData`, `zFaceData`.

#### `zFnType`

Source: [bind_core_types.cpp](core/src/bind_core_types.cpp).

Enum values: `zInvalidFn`, `zPointsFn`, `zMeshFn`, `zGraphFn`, `zMeshFieldFn`, `zPointFieldFn`, `zParticleFn`, `zMeshDynamicsFn`, `zGraphDynamicsFn`, `zComputeMeshFn`, `zNurbsCurveFn`, `zArcFn`, `zPlaneFn`.

#### `zPlanarSolverType`

Source: [bind_core_types.cpp](core/src/bind_core_types.cpp).

Enum values: `zQuadPlanar`, `zVolumePlanar`.

#### `zFnGraph`

Source: [bind_fn_graph.cpp](core/src/bind_fn_graph.cpp).

Properties: `planarGraph`, `graphNormal`.

Methods: `setObject`, `object`, `hasObject`, `getType`, `getBounds`, `clear`, `create`, `createFromMesh`, `numVertices`, `numEdges`, `computeEdgeColorfromVertexColor`, `computeVertexColorfromEdgeColor`, `averageVertices`, `setVertexPositions`, `setVertexColor`, `setVertexColors`, `setEdgeColor`, `setEdgeColors`, `setEdgeWeight`, `setEdgeWeights`, `getVertexPositions`, `getVertexColors`, `getVertexWeights`, `getEdgeColors`, `getEdgeWeights`, `getCenter`, `getCenters`, `getHalfEdgeLengths`, `getEdgeLengths`, `getEdgeData`, `getDuplicate`, `getGraphMesh`, `setTransform`, `setScale`, `setRotation`, `setTranslation`, `setPivot`, `getTransform`.

#### `zFnMesh`

Source: [bind_fn_mesh.cpp](core/src/bind_fn_mesh.cpp).

Methods: `setObject`, `object`, `hasObject`, `getType`, `getBounds`, `clear`, `reserve`, `create`, `numVertices`, `numEdges`, `numPolygons`, `getTopologicalLaplacian`, `checkPointInConvexHull`, `computeEdgeColorfromVertexColor`, `computeVertexColorfromEdgeColor`, `computeFaceColorfromVertexColor`, `computeVertexColorfromFaceColor`, `smoothColors`, `computeVertexNormalfromFaceNormal`, `computeMeshNormals`, `averageVertices`, `makeConvexHull`, `isTriMesh`, `isQuadMesh`, `setVertexPositions`, `setVertexColor`, `setVertexColorsfromScalars`, `setVertexColors`, `setVertexWeight`, `setVertexWeights`, `setFaceColor`, `setFaceColors`, `setFaceColorOcclusion`, `setFaceNormals`, `setEdgeColor`, `setEdgeColors`, `setEdgeWeight`, `setEdgeWeights`, `getVertexPositions`, `getVertexNormals`, `getVertexColors`, `getVertexWeights`, `getEdgeColors`, `getEdgeWeights`, `getFaceNormals`, `getFaceColors`, `getCenter`, `getCenters`, `getDualMesh`, `getDualGraph`, `getRainflowGraph`, `getMeshTriangles`, `getMeshVolume`, `getMeshFaceVolumes`, `getPrincipalCurvatures`, `getGaussianCurvature`, `getPlanarityDeviationPerFace`, `getEdgeDihedralAngles`, `getHalfEdgeLengths`, `getEdgeLengths`, `getVertexAreas`, `getPlanarFaceAreas`, `getPolygonData`, `getMatrices_trimesh`, `getMatrices_quadmesh`, `getEdgeData`, `getDuplicate`, `splitMesh_Mixed`, `splitMesh_Quad`, `getIsoContour`, `getIsoMesh`, `getIsoMesh_mixed`, `getIsobandMesh`, `triangulate`, `subdivide`, `smoothMesh`, `extrudeMesh`, `extrudeVariableMesh`, `extrudeBoundaryEdge`, `setTransform`, `setScale`, `setRotation`, `setTranslation`, `setPivot`, `getTransform`.

#### `zObjGraph`

Source: [bind_obj_graph.cpp](core/src/bind_obj_graph.cpp).

Methods: `__repr__`.

#### `zObjMesh`

Source: [bind_obj_mesh.cpp](core/src/bind_obj_mesh.cpp).

Methods: `__repr__`.

#### `ObjMeshScalarField`

Source: [module.cpp](core/src/module.cpp).

Default construction only; no additional methods are registered.

### zspace.toolset.slicer

#### `zTsSlicing`

Source: [bind_ts3dp.cpp](toolsets/zTs3DP/src/bind_ts3dp.cpp).

Methods: `setInputMesh`, `setBracingGraph`, `setSlotAnchor`, `clearSlotAnchor`, `hasSlotAnchor`, `getSlotAnchor`, `getFaceTags`, `getFaceWidthIds`, `setSlicerMode`, `getSlicerMode`, `computeSlicer`, `computePrintPath`, `computeBasicStructure`, `computePostProcess`, `isLoaded`, `getStatus`, `getSettings`, `getWidthGroupCount`, `getPrintSynSetting`, `get2DPrintPaths`, `get3DPrintPaths`, `getPrintHeightDomain`, `getPrintWidthDomain`, `getPrintLength`, `getPrintVolume`, `hasPrintStatistics`, `getResultantPrintPath`, `getIsoMeshes`, `getIsoMeshes2D`, `getScalarField`, `getPathFields`, `getWidthFields`, `getTrackFields`, `getInputMesh`, `getBracingGraph`.

#### `SlicerMode`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Enum values: `Topology`, `Harmonic`.

#### `SlicerSettings`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Properties: `BoundaryWidth`, `BracingWidth`, `BracingCurveWidths`, `OverlapWidth`, `DoubleBoundary`, `SlotOffset`, `BracingSlot`, `FieldScaleUp`, `FieldCellSize`, `LayerSpacing`, `_HTraceStepFraction`, `_HTraceMaxSteps`, `_HEqualiseIterations`, `_HArapIterations`, `_HUnrollStart`.

#### `VertexTag`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Enum values: `Interior`, `Bottom`, `Top`, `Side`.

#### `PrintSynSetting`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Properties: `ResampleLength`, `FeatureAngleTol`, `SeamOffset`, `SeamSize`.

#### `ParameterGraphVertex`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Properties: `v0`, `v1`, `v2`, `u`, `v`, `w`.

#### `ParameterGraph`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Properties: `vertices`, `edgeConnects`.

Methods: `empty`, `clear`.

#### `PrintPath`

Source: [bind_types.cpp](toolsets/zTs3DP/src/bind_types.cpp).

Properties: `graph`, `widths`, `heights`, `normals`, `trackIds`.

#### `USDdoc` (in `zspace.core`)

Source: [module.cpp](core/src/module.cpp), implemented by Core IO `USDdoc`.

Methods: `supported`, `open`, `close`, `readMetadata`, `readPrim`,
`readAttribute`, `readCustomLayerData`, `readPoints`, `applyEdits`, `save`, `text`.
Queries return native Python dictionaries/lists, not JSON strings. Vector arrays
are flat XYZ. `readPoints` applies world transforms but does not convert units;
`readAttribute` expands indexed primvars. `applyEdits` accepts dictionaries of
typed values and commits only after every change succeeds. The binding converts
these directly to C++ structures; it never calls JSON loads/dumps. Numeric
sequences are copied to typed vectors, so this is not a zero-copy API.
The old query/edit names are removed. Rebuild the Core DLL and Python extension
together. Saving is explicit; concurrent mutation is unsupported.

With USD enabled, `zspace.toolset.slicer` exposes `SlicerUSDInput`,
`importUSD(path)` and `exportUSD(print_path, input_path, output_path, layercount, blockid=0)`.
These call the same Toolsets C++ IO used by Alice and WASM, including units,
settings defaults and output statistics. Keep the solver owning the PrintPath alive during export.

`zspace.io.slicing_usd.read_input(path)` is a thin buffer adapter for Blender.
`write_output(input_path, output_path, result)` exports the native PrintPath retained
in the Blender compute result. There is no separate Python USD schema parser.

## Module constants

`zspace.core` provides color constants `zRED`, `zGREEN`, `zBLUE`, `zWHITE`,
`zGREY`, `zBLACK`, `zMAGENTA`, `zYELLOW`, and `zCYAN`.

## Slicer usage notes

`zTsSlicing` is the exposed facade for Topology and Harmonic modes. `TopologySlicer`, `HarmonicSlicer`, `zPrintSyn` and utility namespaces are not independently bound. `ObjMeshScalarField` is exposed as a container with a default constructor, not a full field editing API. Enum names do not imply corresponding classes are bound.

Supply a mesh, face tags/width IDs, bracing as needed, and an explicit slot anchor. Configure the live settings returned by `getSettings()` and `getPrintSynSetting()`. The compute bindings accept `use_parallel=True` and release the Python GIL while native processing runs. Avoid concurrent mutations of the same engine or its inputs.

Direct geometry and dimensional settings should use the slicer's internal millimetres. Ratios and angles are not length-scaled. USD-enabled input performs its own unit conversion. Inspect the return status and `getStatus()` before consuming results. Input object and returned-reference lifetime behavior is defined in the binding source; keep the owning engine alive while using its results.

For a complete headless pipeline example, see [scripts/smoke_headless.py](scripts/smoke_headless.py). It is not executed automatically by the build script. Source changes require rebuilding before installed packages reflect them.

## Shared implementation sources

Python enables `ZSPACE_TOOLSETS_SLICING_ONLY`. Its native library and
`wasm/slicing` use `cmake/SlicingSources.cmake` to select the same
`src/zToolsets/fabrication/slicing` sources. Harmonic files are included when
TetGen is enabled. Other toolsets are excluded from the Python build.
Core configuration, platform bridges and linking remain platform-specific.
USD is implemented by Core IO on native and WASM targets. Toolsets IO
interprets the slicer schema consistently across native, Python and WASM builds.

USD SDK selection and runtime deployment are owned by Core. Configure Core
`cmake/OpenUSD.local.cmake` or its platform environment override once.

`readPoints` returns `std::vector<std::array<double, 3>>` in C++ and a list
of XYZ triples in Python. The WASM bridge flattens these triples for its numeric
buffer contract, so the browser API still returns a flat coordinate array.


### Blender geometry-only build

With `ZSPACE_PYTHON_WITH_USD=OFF`, Core IO is not built or linked, `USDdoc` is not
bound, and no `zspace.io` adapter is generated in a fresh package. Core/Interface
and Toolsets remain native DLLs with their separate `_core`/`_slicer` PYD bindings.
The Blender extension requires USD enabled; see sibling
`zspace_blender/build_extension.ps1`. Use a fresh build directory when switching
from an older USD-enabled deployment to avoid packaging stale IO files.

## Core standalone and combined builds

Core can be built independently using `zspace_core/python/build_headless.ps1`. This Toolsets entry point builds native Core and Toolsets together, and includes the Core-owned binding target to produce `_core` and `_slicer` in one package. No preinstalled Core SDK is required. `ZSPACE_CORE_SOURCE` selects the matching Core checkout. Headless Python configuration is shared through Core `python/HeadlessPython.cmake`. Use a fresh build directory after this binding relocation.

## NumPy result extraction

`PrintPath.to_numpy()` returns a dict of `points`, `edges`, `normals`, `widths`, `heights`, and `trackIds`. XYZ arrays have shape `(N,3)`; other arrays are one-dimensional. Arrays own independent copies and remain valid after the path/solver is destroyed. Graph order is preserved; traversal ordering is an application responsibility. Core also provides NumPy Mesh/Graph creation and extraction. Requires NumPy at runtime.

## Shared Python runtime directory

The package stores native DLLs only in `zspace/runtime/`. `core/` and
`toolset/slicer/` contain their Python initializers and ABI-tagged PYD modules.
The root initializer registers the runtime DLL directory before importing Core.
USD resources are stored once in `runtime/usd/` and `runtime/plugin/`.
The `zspace/io/` Python adapter remains part of the package.
Windows combined builds deploy the OpenMP runtime automatically; override
`ZSPACE_PYTHON_OPENMP_RUNTIME` when the MSVC redistributable is elsewhere.
