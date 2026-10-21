# zSpace Python bindings

Bindings live in `zspace_toolsets/python`. C++ algorithms remain in the v3
`include/zspace` and `src` trees. Python imports remain `zspace.core` and
`zspace.toolset.slicer`.

## Layout

- `build_headless.ps1`: Windows build entry point for Python 3.10 and 3.13.
- `CMakeLists.txt`: headless core, v3 toolsets and Python package configuration.
- `core/src`: geometry bindings.
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
`-WithUSD -OpenUSDRoot C:/SDK/OpenUSD` enables its bindings through the v3 helper.
For USD builds, deploy the matching SDK runtime DLLs and plugin resources using
the SDK's runtime layout/search paths; these are not automatically packaged.

CMake overrides include `ZSPACE_CORE_SOURCE`, `ZSPACE_TOOLSETS_SOURCE`,
`ZSPACE_PYTHON_BUILD_SLICER`, `ZSPACE_PYTHON_WITH_USD`,
`ZSPACE_TOOLSETS_WITH_TETGEN`, and `ZSPACE_TETGEN_ROOT`.
Core display, interop, tests and OpenGL are disabled. Core IO is built because
v3 toolsets links it; core OpenUSD support remains disabled.

## Package and API compatibility

Set `PYTHONPATH` to the generated `package` directory. `_core.pyd` links native
core/interface DLLs; `_slicer.pyd` links `zSpace_Toolsets.dll`, with `zSpace_IO.dll`.
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

#### `SlicerUSDInput` (USD builds only)

Source: [module.cpp](toolsets/zTs3DP/src/module.cpp).

Properties: `block`, `face_tags`, `face_width_ids`, `bracing_curves`, `bracing_width_ids`, `bracing_print_widths`, `slotpoint`, `settings_json`, `meters_per_unit`, `up_axis`.

Methods: `read_input_usd`, `export_printpath_usd`.

## Module constants and optional functions

`zspace.core` provides color constants `zRED`, `zGREEN`, `zBLUE`, `zWHITE`, `zGREY`, `zBLACK`, `zMAGENTA`, `zYELLOW`, and `zCYAN`. Vector multiplication supports a vector dot product, scalar multiplication and matrix multiplication; `^` is the cross product.

`zspace.toolset.slicer.HAS_USD` reports the compile-time USD capability. When enabled, the module additionally exposes:

- `read_input_usd(file) -> SlicerUSDInput`.
- `export_printpath_usd(path, input_file, output_file, layercount, blockid=0)`.

Both helpers raise an exception on reported failure. They are not registered when USD is disabled.

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
The separate WASM USD reader intentionally compiles only its IO utility subset.
