# C++ slicer and PrintSyn

`zSpace::zTsSlicing` is the maintained toolset API for topology and harmonic
slicing. It uses the geometry kernel from `zspace_core` and supplies results to
the desktop and browser applications. This guide covers slicer only.

Start with the [project overview](https://github.com/GitZHACODE/zspace_core/blob/feature/zspace-v3-gitbook/docs/gitbook/slicer/README.md),
then read the [input contract](input.md) and [build configurations](build.md).
Cross-repository mismatches are tracked in the
[compatibility register](https://github.com/GitZHACODE/zspace_core/blob/feature/zspace-v3-gitbook/docs/gitbook/slicer/compatibility.md).

## Pipeline

Topology slicing uses a tagged all-quad mesh with a connected floor patch and
regular columns from floor to top. Harmonic slicing uses a tetrahedralized volume
and requires TetGen support. Query `supportsSlicerMode()` before selecting
harmonic mode in a build where it may be disabled.

Both modes construct layer geometry, charts and scalar fields. Path extraction
maps field contours back to 3D. PrintSyn then handles path resampling, seam
processing, projection and linking.

| Method | Result boundary |
| --- | --- |
| `computeBasicStructure()` | Input validation and initial topology or tetrahedral structure |
| `computeSlicer()` | Layer geometry, charts and fields |
| `computePrintPath()` | Extracted and mapped paths without PrintSyn |
| `computePostProcess()` | Processed per-layer paths and the combined path |

## C++ integration

Include `<zspace/slicing.h>` and link against `zspace::toolsets`. Given a valid
mesh, per-face tags and width IDs, optional bracing, and a slot anchor:

```cpp
#include <zspace/slicing.h>
#include <stdexcept>

zSpace::zTsSlicing session;
session.setInputMesh(mesh, faceTags, faceWidthIds);
session.setBracingGraph(bracing); // omit when no bracing is needed
session.setSlotAnchor(anchor);
session.getSettings().LayerSpacing = 1.0f;
session.getSettings().BoundaryWidth = {0.5f};
session.getSettings().FieldCellSize = 0.05f;
session.getPrintSynSetting().ResampleLength = 0.1f;

if (!session.computeSlicer(false) || !session.computePostProcess(false))
    throw std::runtime_error(session.getStatus());

const auto& layers = session.get3DPrintPaths();
const auto& combined = session.getResultantPrintPath();
```

The default mode is topology. `false` requests sequential execution. Input
variables in this fragment must be supplied by the application; it is not a
standalone program or a benchmark.

The session owns copies of its mesh, tags, width IDs, bracing and anchor. Direct
`TopologySlicer` / `HarmonicSlicer` operators borrow their inputs. Calling
`setSlicerMode()` recreates the operator and invalidates references to its results,
even when the requested mode is unchanged.

## Outputs

Each print path contains graph geometry, widths, heights, normals and track IDs.
Use `get3DPrintPaths()` for layer paths and `getResultantPrintPath()` for the
combined path. Check return values and status before consuming results; partial
diagnostic geometry is not a completed print result.

Harmonic heights use the local normal-spacing estimate `deltaW / |grad(w)|`.
They are not exact distances between curved layers or a maximum guaranteed by
`LayerSpacing`. `zTsUtils::buildPrintMesh()` lofts rectangular sections using
path widths, heights and normals. Native USD export can preserve the source
input alongside print-path and print-mesh output.

## Sources

Source baseline: `4a8329c72415917f871f1ff1bc825f1d39155fbf`.

- [Slicing API](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/include/zspace/zToolsets/fabrication/slicing/zTsSlicing.h)
- [Tagged slicing implementation guide](https://github.com/GitZHACODE/zspace_toolsets/blob/4a8329c72415917f871f1ff1bc825f1d39155fbf/docs/architecture/tagged-slicing.md)
