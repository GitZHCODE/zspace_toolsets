# Carbcomn Grasshopper JSON Exporter

Use `CarbcomnJsonExporter_GhPython.py` in a Grasshopper GhPython component to write `blockMesh_<id>.json` files for the current Carbcomn slicer.

The current Carbcomn path is simplified to non-planar regular wall blocks. The slicer no longer uses `BlockType`; do not write it.

## Required Data

The JSON has two parts:

- A zSpace half-edge mesh, used by `zFnMesh::from(path, zJSON)`.
- Carbcomn slicing metadata, read by `zTsCarbcomn::readJSON`.

## zSpace Mesh Fields

`Vertices`

Stores one outgoing half-edge id per vertex. Carbcomn uses the reconstructed zSpace mesh topology for all later walking operations. This is not just display geometry; invalid half-edge connectivity will break medial walks, feature-edge coloring, geodesic contours, and unrolling.

`Halfedges`

Stores each half-edge as `[previousHalfedgeId, nextHalfedgeId, endVertexId, faceId]`. Half-edges are expected in paired order, so half-edge `0` pairs with `1`, `2` with `3`, and so on. This is the most important topology field.

`Faces`

Stores one half-edge id per face. This is used to reconstruct face loops, normals, and downstream mesh traversal.

`VertexAttributes`

Stores per-vertex data as `[x, y, z, nx, ny, nz, r, g, b]`. Positions define the input block geometry. Normals and colors are currently secondary, but zSpace expects them in the JSON.

`HalfedgeAttributes`

Stores edge colors as `[r, g, b]` per half-edge. Carbcomn recolors feature/corner edges after loading, so exported values can be black.

`FaceAttributes`

Stores per-face data as `[nx, ny, nz, r, g, b]`. Face normals are useful for downstream checks and exports. Face colors are secondary.

## Carbcomn Metadata

`IsPlanar`

Must be `false`. Carbcomn currently forces non-planar slicing, but the key is kept for compatibility with older JSON files.

`IsCorner`

Must be `false`. The current cleaned path assumes regular wall topology, not corner/pentagon topology.

`MedialStartEnd`

`[startVertexId, endVertexId]`. These vertex ids define the medial/spine direction on the input mesh. Carbcomn uses them to:

- find the start half-edge with `util_getStartHalfEdge`,
- build `o_MedialGraph`,
- seed the regular slice mesh coloring,
- compute vertical loops for non-planar sectioning.

`FeaturedNumStrides`

Odd-length list of positive integer stride counts. This describes how many mesh edge steps exist between successive feature/corner bands around the section topology. Carbcomn uses it in `compute_SliceMesh_Regular` to color inner/outer corners and feature edges, then later uses those colors to build trims, bracing, offsets, and SDF fields.

The list should have a clear middle entry. The middle stride is treated as the bridge between inner and outer sides.

`StartCornerVID`

Mesh vertex id for the inner start corner. Carbcomn uses this vertex to locate the start corner edge and assign the `_col_in_corner_st` / `_col_out_corner_st` topology colors during slice-mesh setup.

## GhPython Inputs

Create these inputs on the component:

- `M`: Rhino mesh
- `folder`: output folder path
- `block_id`: integer id
- `medial_start`: mesh vertex id
- `medial_end`: mesh vertex id
- `feature_strides`: list of positive integers
- `start_corner_vid`: mesh vertex id
- `write`: boolean trigger

Outputs:

- `path`
- `message`

## Slicing Workflow Conclusion

The current Carbcomn workflow should be understood as a topology-driven non-planar wall slicer. The input mesh is not treated as an anonymous display mesh. It is expected to contain stable vertex ids, valid half-edge connectivity, and a regular wall section layout that the slicer can walk consistently from the medial start vertex to the medial end vertex.

The workflow starts in Rhino or Grasshopper with a raw block mesh. Before export, the user must identify three pieces of slicing metadata: the medial start/end vertices, the ordered feature stride counts, and the start corner vertex. These values define how Carbcomn interprets the mesh topology. No block type, corner flag, or left/right plane metadata is required for the current cleaned non-planar path.

When the JSON is loaded, `zTsCarbcomn::readJSON` reconstructs the zSpace mesh from the half-edge data. The mesh topology is then used to find the start half-edge between `MedialStartEnd[0]` and `MedialStartEnd[1]`. From this, Carbcomn builds the medial graph, which acts as the spine of the block. This spine controls the longitudinal direction of the slicing workflow.

After the medial graph is built, Carbcomn creates the regular slice mesh with `compute_SliceMesh_Regular`. `FeaturedNumStrides` tells the slicer how far to walk between important section bands. During this stage, edges and vertices are colored into semantic groups such as inner corner, outer corner, inner feature, outer feature, and start corner bands. These colors are not cosmetic; later stages use them as topology labels.

The non-planar sectioning stage then uses the medial graph and colored slice mesh to generate section graphs along the block. The current path uses geodesic / mesh-based contour logic rather than planar start/end planes. This is why `LeftPlanes` and `RightPlanes` are no longer part of the required JSON schema.

Once section graphs exist, Carbcomn computes print blocks and section meshes. It evaluates layer heights, checks section geometry, and builds trim-related graph features. The same topology colors assigned earlier are used to decide where hard trim features, soft trim features, bracing, slots, seam alignment, and interior split features should be generated.

The SDF stage builds scalar fields from the non-planar section data and trim graphs. These fields are used to produce fabrication geometry and validate spacing or layer-height constraints. Because the scalar-field construction depends on earlier section and trim graph placement, incorrect topology metadata at export time will usually appear later as failed section walks, missing trims, bad offsets, or invalid SDF checks.

The final export stage writes the generated geometry for visualization and downstream use. In the current USD export, related elements are grouped into higher-level prims such as section meshes, trim bracing, soft trim features, hard trim features, contours, and other repeated element groups. This makes the result easier to inspect because each generated feature family can be isolated or hidden as a group.

In short, the required JSON is minimal but strict. The mesh provides the geometry and walkable half-edge topology. `MedialStartEnd` provides the slicing direction. `FeaturedNumStrides` provides the section topology pattern. `StartCornerVID` anchors the feature-coloring logic. Everything after loading depends on those four pieces being consistent with the final exported mesh.

## OBJ Feasibility

Using OBJ as the raw geometry source is feasible, but not as a direct drop-in for the current Carbcomn JSON.

OBJ can store vertices and faces, but it does not store the Carbcomn metadata: `MedialStartEnd`, `FeaturedNumStrides`, and `StartCornerVID`. OBJ also does not preserve zSpace half-edge ids, and the slicer depends on stable vertex ids and edge-walk topology.

Recommended approach:

- Use OBJ only as the raw mesh input.
- Import the OBJ into Rhino or Grasshopper.
- Assign or compute the required Carbcomn metadata in Grasshopper.
- Export the final `blockMesh_<id>.json` with this GhPython script.

Direct OBJ loading inside Carbcomn would require code changes:

- call `zFnMesh::fromOBJ` or equivalent before slicing,
- add a sidecar metadata file or JSON wrapper for the Carbcomn-specific fields,
- verify that OBJ vertex ordering matches the ids used by `MedialStartEnd` and `StartCornerVID`,
- rebuild or validate half-edge topology before calling the current slice code.

So OBJ is useful as an interchange format, but the current slicer still needs the Carbcomn JSON metadata layer.

## Notes

- The script triangulates a duplicate of the mesh before writing.
- The mesh must be manifold enough for paired half-edge construction. Duplicate oriented edges or non-manifold edges are rejected.
- `MedialStartEnd`, `FeaturedNumStrides`, and `StartCornerVID` must refer to the final mesh vertex ids used by the exported mesh.
