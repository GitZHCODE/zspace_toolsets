# Isocontour Extraction: Toolsets and Core

This comparison describes the current local implementations inspected on 2026-10-23. It compares source behavior, not measured performance.

## Functions and call paths

Toolsets:

```text
extractPrintPaths()
  -> extractRawIsoContour()
  -> extractTrackPaths()
  -> validate paths and select the longest usable closed loop
```

Core:

```text
zFnMeshScalarField::getIsocontour()
  -> read field values
  -> convert face values to vertex values when necessary
  -> zFnMesh::getIsoContour()
  -> construct the output zObjectGraph
```

Source references:

- [Toolsets Field.cpp](../src/zToolsets/fabrication/slicing/utils/Field.cpp): `extractRawIsoContour()` and `extractPrintPaths()`.
- [Core zFnMeshField.cpp](../../zspace_core/src/zInterface/functionsets/zFnMeshField.cpp): `zFnMeshFieldBase<zScalar>::getIsocontour()`.
- [Core zFnMesh.cpp](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp): `zFnMesh::getIsoContour()` and `interpolateCorner()`.

Core links assume that `zspace_core` and `zspace_toolsets` are sibling repositories.

## Behavior comparison

| Aspect | Toolsets `extractRawIsoContour()` | Core `getIsocontour()` |
| --- | --- | --- |
| Scalar input | Explicit `pathVals` array | Values stored in the field object |
| Value domain | Per vertex only | Per vertex, or per face converted to vertices by averaging adjacent face values |
| Isovalue | Fixed at zero | Configurable `inThreshold` |
| Polygon processing | Fan-triangulates each face, then processes triangle sign cases | Finds intersections on polygon boundary edges and pairs them in traversal order |
| Shared contour vertices | Identifies crossings by source edge IDs, or source vertex IDs for exact endpoint crossings | Uses a coordinate-based position map controlled by `precision` |
| Short segments | Rejects a segment when both endpoints resolve to the same ID; no length-based rejection in this function | Discards segments whose length is at most `distTolerance` |
| Width data | Linearly interpolates `widthVals` at every generated contour vertex | No print-width output |
| Track data | Interpolates track values and classifies them into 0 or 1 | No track-ID output |
| Color data | No color interpolation | The Mesh function interpolates colors, but the Field wrapper does not assign the returned colors to its graph |
| Output | Positions, edge connectivity, widths and track IDs | A constructed `zObjectGraph` |
| Invalid input size | Clears outputs and returns if path/width sizes do not match the vertex count | Throws when field values match neither vertices nor faces |
| Cleanup | Deferred to subsequent Toolsets functions | Position merging and short-segment rejection occur during extraction; no guarantee of a single closed loop |

The Core Field wrapper passes `precision` to the Mesh extraction, but constructs the final graph using `PRECISION`.

## Interpolation and zero handling

Toolsets extracts the zero crossing using:

```cpp
t = abs(s1 - s0) > 1e-9
    ? clamp(-s0 / (s1 - s0), 0.0, 1.0)
    : 0.5;
```

It uses the same `t` for position and width interpolation. If track data has the expected vertex count, it interpolates track values and outputs `1` when the result is at least `0.5`, otherwise `0`. Missing or incorrectly sized track data produces zero track values in this raw function; the enclosing `extractPrintPaths()` performs stricter size validation.

Triangle classification uses `scalar > 0`, so zero belongs to the nonpositive side. Crossings with exactly `t == 0` or `t == 1` use the source vertex ID as their sharing key.

Core uses `scalar - threshold`, explicitly includes a corner when its absolute difference is at most `1e-12`, and interpolates an edge when the endpoint differences have opposite signs. Its interpolation helper uses a `1e-12` denominator threshold with a zero fallback parameter.

These rules differ for exact-zero vertices, nearly equal values and other degenerate cases. The scalar thresholds should not be confused with geometric distance tolerances.

## Polygon topology

For a quad, Toolsets introduces a diagonal from the first vertex and extracts the contour on the resulting two triangles. A crossing on that diagonal can introduce an additional contour vertex.

Core processes only the polygon perimeter. It gathers boundary intersections and connects consecutive pairs. For a face with alternating positive and negative corner values, its connections can differ from the fan-triangulated result.

Consequently, the two functions can produce different vertex counts, segment geometry or connectivity from the same polygon mesh and scalar values. Matching inputs and thresholds alone does not guarantee matching outputs.

## Vertex sharing and small features

Toolsets shares crossings through source topology. Two geometrically coincident points on unrelated source edges are not automatically welded by `extractRawIsoContour()`.

Core shares positions through coordinate precision, independently of source edge identity. It also rejects short segments before adding their endpoints to the output position map. Position merging and segment rejection therefore use different criteria.

In the Toolsets pipeline, `extractPrintPaths()` subsequently computes a weld tolerance from the field grid and passes it to `extractTrackPaths()`. Raw extraction and later cleanup must be distinguished when comparing complete results.

Neither raw extraction nor Core graph construction guarantees that the result is a single closed, unbranched loop.

## Suitability for the current slicer

The Toolsets implementation directly preserves the correspondence between contour vertices, print widths and track IDs. It also provides the triangulation and topology-sharing behavior expected by the subsequent path cleanup.

Replacing it with the existing Core function requires decisions and validation for:

1. Width interpolation and binary track classification at generated vertices.
2. Polygon triangulation and ambiguous multi-crossing faces.
3. Exact-zero and near-zero scalar handling.
4. Position sharing, short segments and weld tolerances.
5. Downstream loop assembly and longest-loop selection.

For current behavior, retain `extractRawIsoContour()` unless these differences are intentionally reconciled. A future shared Core extractor could expose source-edge IDs and interpolation parameters so callers can interpolate arbitrary attributes without making Core depend on slicer-specific width or track concepts.

## Performance limits of this comparison

No benchmark was performed. Toolsets copies grid positions and polygon arrays before extraction and uses integer edge keys. Core reads mesh storage and uses a string-keyed coordinate map, then constructs a graph. These implementation differences alone do not establish which complete pipeline is faster; equivalent output and cleanup requirements must be included in a fair benchmark.
