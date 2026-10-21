# zspace_toolsets

Modern zSpace toolsets are built as focused libraries on top of `zspace_core`.

`zTsForcePolyhedra` constructs disjointed force cells from graphs using EGI and
area pursuit, then couples the form/force diagrams through weighted geometric
3DGS reciprocation. Supports, nodal self weight and branching graphs are supported.
See [Force Polyhedra](docs/architecture/force-polyhedra.md).

`zTsParametricStairs` constructs variable-width graph ribbons with automatic
end supports and compression-only, fixed-XY Best Fit FDM. The Alice sketch
updates automatically as its graph and width controls change.
See [Parametric Stairs](docs/architecture/parametric-stairs.md).

The pre-v3 header/source/project layout has been moved to `legacy/pre-v3`. New work should use:

```text
include/zspace/zToolsets/
src/zToolsets/
examples/
tests/
```

The maintained slicer is `zTsSlicing` under `fabrication/slicing`, providing tagged Topology/Harmonic slicing, field generation, contour extraction and PrintSyn. Include `<zspace/slicing.h>` or `<zspace/toolsets.h>`.

`zTsVault` provides Best Fit Force Density modelling with support vertices,
tension/compression assignments, vertex masses and bounded force densities.
See [Vault API, build and viewer workflow](docs/architecture/vault-best-fit-fdm.md).

`zTsSlicing` adds the desktop tagged topology/harmonic slicing and PrintSyn
pipeline through `<zspace/slicing.h>`. See [tagged slicing](docs/architecture/tagged-slicing.md)
for input tags, optional TetGen/OpenUSD/OpenMP dependencies, and native tests.

Build:

```bat
scripts\build_toolsets.bat
```

Browser transfer APIs and compatibility are documented in
[typed WASM transfers](docs/architecture/typed-wasm-transfers.md).
