# Input geometry, tags and units

The compute API accepts geometry and settings; it performs no file IO. Supply a
mesh, one tag per polygon, a slot anchor, and optional face-width IDs and bracing
in the same coordinate system. Coordinates and length settings must use the
same unit. The native USD adapter normalizes them to millimetres.

## Face tags and width groups

| Tag | Value | Meaning |
| --- | --- | --- |
| Bottom | 1 | Floor patch |
| Top | 2 | Top patch |
| Side | 3 | Side faces |

Despite the `VertexTag` type name in the source, the input table contains one
entry per face, in polygon order. Side-face width IDs index `BoundaryWidth`.
The native tagged USD reader requires `-1` for top/bottom width IDs and validates
side IDs against the supplied width table.

Topology mode requires all-quad faces and regular floor-to-top columns. Harmonic
mode requires a volume suitable for tetrahedralization and a TetGen-enabled build.
Both operators require the slot anchor. Bracing supplied directly to the compute
session is optional.

## Native tagged USD schema

The current native adapter `zTsIO::readSlicerInputUSD` expects these default paths:

| Prim | Content |
| --- | --- |
| `/AG_block` | Polygon Mesh |
| `/Slicer_bracing` | Linear BasisCurves, periodic or nonperiodic |
| `/Slicer_slot` | Slot point geometry |
| `/Slicer_data` | Scope containing `block:FaceTag`, `block:widthid` and `bracing:widthid` tables |
| `/Slicer_setting` | Scope containing `slicer:*` and `printSyn:*` settings attributes |

The geometry paths can be supplied explicitly to the reader. Tag and setting
Scopes use the names above. Each bracing width ID selects a width for its curve.
The reader validates geometry, table lengths, IDs and settings before returning
a complete input. Mesh hole faces are rejected.

Object transforms are applied. Geometry and dimensional settings are multiplied
by `metersPerUnit * 1000`; ratios, angles, tags and IDs are unchanged. The original
`MetersPerUnit` is retained for export. `UpAxis` is recorded; this reader does
not promise the automatic Y-up rotation described by the older browser guide.

A checked-in [Alice 2 input fixture](https://github.com/GitZHACODE/zspace_alice2/blob/d8d3bd4a1e03cd19ecd16e4231b9054efde79184/alice2/verification/fixtures/slicing-block.usda)
shows the complete current schema. Earlier descriptions using `/block`,
`/bracing`, `/slotpoint` and root-layer setting dictionaries describe a different
input contract. Check the deployed importer before transferring a file between
desktop and browser versions.

## Print export

`zTsIO::exportPrintPathUSD()` copies the source stage and replaces `/TDP_gcode`,
`/TDP_printmesh` and `/TDP_printpath`. Internal millimetric geometry, widths and
heights are converted back to source stage units. Pass processed per-layer paths
through `layerPaths` to keep the visualization separate across layers and avoid
lofting inter-layer travel connectors.

## Sources

- [Native USD reader and writer](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/src/zToolsets/fabrication/slicing/utils/IO.cpp)
- [Input and settings declarations](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/include/zspace/zToolsets/fabrication/slicing/utils/Utils.h)
