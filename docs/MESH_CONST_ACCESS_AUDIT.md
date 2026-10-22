# Mesh Read-Access and Cache Invalidation Audit

Date: 2026-10-21. Scope: current local zspace_core v3 sources and TopologySlicer callers. This is a source audit, not a timing benchmark. No C++ changes or builds were performed for this report. Links are relative to this document and assume sibling core/toolsets repositories. Line numbers are a snapshot.

## 1. Cause before proposed changes

The mesh object keeps a face-list representation and a lazily built half-edge topology. The non-const storage accessor treats access as potential mutation. Several query implementations accidentally use this accessor because their internal mesh pointer is mutable.

Relevant implementation: [zMeshObjectStorage.h](../../zspace_core/src/zInterface/objects/zMeshObjectStorage.h#L73).

| Accessor | Preparation | Sets a new dirty flag? |
|---|---|---|
| `get(zObjectMesh&)` | `buildTopology()`; synchronizes an already dirty face-list | Yes: `faceListDirty = true` |
| `get(const zObjectMesh&)` | `buildTopology()`; synchronizes an already dirty face-list | No |
| `read(const zObjectMesh&)` | `syncFaceList()` only when needed | No |
| `edit(zObjectMesh&)` | Synchronizes face-list, then discards topology | Yes: topology and position lookup invalidation |

`buildTopology()` checks `faceListDirty` BEFORE returning an existing clean topology. `syncFaceList()` enumerates active vertices and faces, copies geometry/attributes, calls `rebuilt.set()`, rebuilds face-list edges, constructs an edge mapping, and transfers edge attributes. This is whole-mesh work, not a constant-time flag update. See [syncFaceList](../../zspace_core/src/zInterface/objects/zMeshObjectStorage.h#L120) and [face-list initialization](../../zspace_core/src/zCore/geometry/detail/zMeshFaceListStorage.cpp#L45).

Distinguish rebuilding the face-list from rebuilding the half-edge topology: repeated non-const reads primarily trigger face-list synchronization; `edit()` can additionally cause a later topology reconstruction. Neither should happen repeatedly during an unchanged read-only traversal.

### Const on the result is insufficient

```cpp
// mesh is mutable: overload resolution chooses the mutable accessor first.
const auto& topology = zMeshObjectStorage::get(mesh);

// A const argument chooses the read-only topology overload.
const auto& topology = zMeshObjectStorage::get(
    static_cast<const zObjectMesh&>(mesh));
```

A trailing `const` on an iterator method also does not make the object behind its mutable pointer const. The storage call itself must use the read-only overload. Adding const public overloads is an API improvement, but is not a prerequisite for correcting existing query implementations.

### Concrete amplification chain

`zItMeshHalfEdge::getNext()` calls `isActive()`, which calls mutable `get()` and marks the mesh dirty. `getNext()` then calls mutable `get()` again, causing a complete face-list synchronization and marking it dirty again. Constructing the returned indexed half-edge iterator also uses mutable `get()`. Even a single next-edge query can therefore cause multiple synchronizations.

If Q queries each trigger whole-mesh work, traversal cost grows with Q times the mesh synchronization cost (including allocations and edge maps), rather than just Q constant-time adjacency lookups. This is a structural cost argument; no measured slowdown factor is claimed.

## 2. Direct query call sites requiring read-only access

The following inventory covers active `zMeshObjectStorage::get()` calls in mesh iterator/functionset implementations. Commented-out code is excluded. Each overload is listed separately. Most fixes can retain the public signature and change only the internal access path.

| Function | Definition | Direct get() lines | Required treatment |
|---|---|---|---|
| `zItMeshVertex::getConnectedHalfEdges` | [zItMesh.cpp:84](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L84) | 86 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshVertex::getConnectedHalfEdges` | [zItMesh.cpp:103](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L103) | 105 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshVertex::getArea` | [zItMesh.cpp:503](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L503) | 539, 542 | Read geometry through const storage; use a local zUtilsCore for getTriangleArea (utility method is not const). |
| `zItMeshVertex::getHalfEdge` | [zItMesh.cpp:560](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L560) | 563 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshEdge::getDihedralAngle` | [zItMesh.cpp:815](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L815) | 825, 826 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshEdge::getHalfEdge` | [zItMesh.cpp:853](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L853) | 857 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshFace::getHalfEdges` | [zItMesh.cpp:991](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L991) | 1009 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshFace::getHalfEdge` | [zItMesh.cpp:1553](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1553) | 1556 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::zItMeshHalfEdge` | [zItMesh.cpp:1816](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1816) | 1821 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::end` | [zItMesh.cpp:1841](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1841) | 1843 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::size` | [zItMesh.cpp:1851](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1851) | 1854 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getVertexPositions` | [zItMesh.cpp:1886](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1886) | 1894 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::onBoundary` | [zItMesh.cpp:1952](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1952) | 1955 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getCenter` | [zItMesh.cpp:1958](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1958) | 1963, 1963 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getVector` | [zItMesh.cpp:1966](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1966) | 1972, 1972 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getSym` | [zItMesh.cpp:1990](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1990) | 1993 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getNext` | [zItMesh.cpp:1998](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1998) | 2001 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getPrev` | [zItMesh.cpp:2006](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2006) | 2009 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getVertex` | [zItMesh.cpp:2014](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2014) | 2017 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getFace` | [zItMesh.cpp:2022](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2022) | 2025 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getEdge` | [zItMesh.cpp:2030](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2030) | 2033 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::getColor` | [zItMesh.cpp:2038](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2038) | 2041 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zItMeshHalfEdge::isActive` | [zItMesh.cpp:2132](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2132) | 2135 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zFnMesh::numHalfEdges` | [zFnMesh.cpp:858](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L858) | 860 | Use const topology access; for position/normal values, read() is also suitable when index correspondence is preserved. |
| `zFnMesh::halfEdgeExists` | [zFnMesh.cpp:878](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L878) | 880 | Query-only intent, but underlying halfEdgeExists is not const-qualified; add a const lookup implementation/overload after auditing its body. Do not use const_cast as a workaround. |
| `zFnMesh::getDualMesh` | [zFnMesh.cpp:1871](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L1871) | 1926, 1955, 1985, 1991, 1996, 1996, 2000, 2002, 2003, 2003, 2010, 2024 | Mixed operation: input mesh reads must be const; output mesh creation/rotation must remain mutable. Use local zUtilsCore; batch output edits and finish them before readback. |
| `zFnMesh::getDualGraph` | [zFnMesh.cpp:2035](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L2035) | 2073, 2090, 2102 | Read source mesh topology/normals through const storage; creating the output graph remains a write to a different object. |
| `zFnMesh::getVertexAreas` | [zFnMesh.cpp:2776](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L2776) | 2804, 2809, 2812 | Read positions through read()/const storage; use local zUtilsCore. This writes output arrays, not the input mesh. |

## 3. Indirectly affected queries

These need call-chain verification after the direct fixes; they do not all need new storage calls. Changing only getNext() while leaving isActive() or indexed constructors mutable is insufficient.

| Class | Query families to verify | Dependency |
|---|---|---|
| `zItMeshVertex` | `getConnectedEdges`, `getConnectedVertices`, `getConnectedFaces` (both overload families), `onBoundary`, `getValence`, neighborhood/curvature/area queries | Connected half-edge traversal and vertex getHalfEdge |
| `zItMeshEdge` | `getVertices`, `getVertexPositions`, `onBoundary`, `getCenter`, `getVector`, `getLength` | Edge getHalfEdge and half-edge navigation |
| `zItMeshFace` | `getHalfEdges` overloads, `getConnectedFaces`, `onBoundary`, traversal-based geometric queries | Face getHalfEdge and half-edge navigation |
| `zItMeshHalfEdge` | `getStartVertex`, `getVertices` overloads, `getConnectedHalfEdges` overloads, `getLength`, `getId`, equality/inequality | isActive, getSym/getVertex, getVector; even comparisons can query active state |
| `zFnMesh` | Dual construction, area queries and iterator-based mesh algorithms | The same primitive iterator reads; audit input reads separately from output writes |

Pure face-list queries such as vertex getPosition(), face getCenter()/getVertices(zIntArray&), and functionset numVertices()/numEdges()/numPolygons(), getVertexPositions() and getPolygonData() already use read(). They can pay for an existing dirty flag, but do not themselves create a new one. Do not replace these with a topology accessor unnecessarily.

## 4. Accesses that must remain writable

These are the remaining active direct get() callers. A blanket replacement would break mutation or pointer semantics.

| Function | Definition | Direct get() lines | Reason |
|---|---|---|---|
| `zItMeshVertex::setHalfEdge` | [zItMesh.cpp:619](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L619) | 626, 628 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshEdge::setHalfEdge` | [zItMesh.cpp:881](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L881) | 887, 889, 890 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshFace::setHalfEdge` | [zItMesh.cpp:1714](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1714) | 1721, 1722 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::deactivate` | [zItMesh.cpp:1857](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L1857) | 1860 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::getRawColor` | [zItMesh.cpp:2044](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2044) | 2047 | Returns a writable pointer; preserve mutation/invalidation semantics. |
| `zItMeshHalfEdge::setId` | [zItMesh.cpp:2052](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2052) | 2055 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::setSym` | [zItMesh.cpp:2059](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2059) | 2061 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::setNext` | [zItMesh.cpp:2068](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2068) | 2075 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::setPrev` | [zItMesh.cpp:2084](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2084) | 2091 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::setVertex` | [zItMesh.cpp:2100](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2100) | 2104 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::setEdge` | [zItMesh.cpp:2110](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2110) | 2114 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zItMeshHalfEdge::setFace` | [zItMesh.cpp:2120](../../zspace_core/src/zInterface/iterators/zItMesh.cpp#L2120) | 2124 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::addVertex` | [zFnMesh.cpp:734](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L734) | 744 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::addEdges` | [zFnMesh.cpp:750](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L750) | 766 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::addPolygon` | [zFnMesh.cpp:775](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L775) | 782 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::addPolygon` | [zFnMesh.cpp:805](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L805) | 807 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::updatePolygon` | [zFnMesh.cpp:813](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L813) | 820 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::addToHalfEdgesMap` | [zFnMesh.cpp:4064](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L4064) | 4066 | Actual mutation; retain writable access. Batch repeated writes where safe. |
| `zFnMesh::removeFromHalfEdgesMap` | [zFnMesh.cpp:4069](../../zspace_core/src/zInterface/functionsets/zFnMesh.cpp#L4069) | 4071 | Actual mutation; retain writable access. Batch repeated writes where safe. |

Other writable APIs use `edit()` rather than `get()`: iterator getRawPosition/getRawNormal and setter families; functionset getRawVertexPositions/getRawVertexNormals/getRawVertexColors/getRawEdgeColors/getRawFaceNormals/getRawFaceColors, and geometry/attribute setters. Their names may contain "get", but the returned memory is writable. Preserve existing semantics; offer separate const-view/copy APIs if needed. A retained raw pointer has lifetime/invalidation constraints and must not be used concurrently with rebuilds.

For mixed functions such as getDualMesh(), use distinct input-read and output-edit scopes. Do not retain references into face-list buffers across a later edit()/syncFaceList(), since these operations can replace storage. Do not retain topology references across edit(), which resets topology. This warning applies when batching accesses too.

## 5. Why TopologySlicer is exposed

[Topology Structure.cpp](../src/zToolsets/fabrication/slicing/zSlicer/topology/Structure.cpp) repeatedly uses the affected queries:

| Helper | Access pattern |
|---|---|
| `_findCapStart` | Vertex incidence, opposite edges, faces and edge lengths |
| `_buildCapMesh` | Face BFS using next/prev/opposite half-edges |
| `_alignBoundaryVertices` | Boundary traversal and indexed vertex access |
| `_traceColumn` | V-loop traversal, edge lengths and vertex valence |
| `_faceTowardColumn` | Adjacent edges/faces and width-group selection |

Harmonic also has some mesh queries, so it is not immune. Its TetGen/array/Eigen-heavy stages make fewer sustained half-edge traversals than this Topology pipeline. Normal Harmonic timing does not disprove the cache problem. ParallelErrors is not on this serial basic-structure path and is not an explanation for this specific repeated synchronization.

Bracing parameterization should be profiled separately: its triangle search and candidate processing are a separate possible cost, not evidence of cache invalidation. The previously removed duplicate sketch computeBasicStructure() invocation was another independent source of overhead.

## 6. Recommended implementation order

1. Fix half-edge isActive(), indexed construction, end()/size(), navigation and value getters together.
2. Fix vertex/edge/face primitive topology queries and connected-half-edge overloads.
3. Replace uses of mesh-owned non-const utility methods with a local zUtilsCore where appropriate. Merely const-casting the mesh argument may not compile for these methods.
4. Resolve the non-const underlying halfEdgeExists API with a verified read-only implementation; inspect [zHalfEdgeTopology.h](../../zspace_core/src/zCore/geometry/detail/zHalfEdgeTopology.h#L69). Do not remove dirty marking globally to accommodate it.
5. Separate source reads and destination writes in dual operations; batch genuine writes without invalidating live references.
6. Rebuild dependent native/WASM artifacts when implementation is authorized. Existing binaries will not change from source edits alone.

The native and WASM consumers share the source defect; it is not specific to OpenMP. Const access does not guarantee thread safety: initialize lazy topology/cache state before parallel readers and prohibit concurrent mutation. A warm, unchanged cache must remain clean throughout a read-only traversal.

## 7. Verification plan (not executed)

- Instrument syncFaceList/buildTopology call counts. Warm the cache, run repeated read-only queries, require zero additional synchronizations/rebuilds.
- Test cold-cache access: allow initial topology creation, then require cache reuse.
- After a real write, verify the face-list reflects the modification and subsequent read-only queries stop synchronizing repeatedly.
- Preserve boundary/opposite-edge semantics, iterator validation, raw-pointer behavior and dual output geometry.
- Compare Topology U-mesh order, V-loop order/lengths, width IDs, bracing parameters and final paths before/after on the same input.
- Time individual basic-structure helpers with the same Release build, input and core revision. Report measurements separately from the static findings above.

## 8. Audit boundaries

A recursive search of current core `src/` found direct zMeshObjectStorage::get() use only in zItMesh.cpp and zFnMesh.cpp. This report inventories those active call sites, reviews the storage implementation, and follows the Topology basic-structure callers. It does not claim every unrelated mesh algorithm, external plugin, or retained raw-pointer user is fully validated. All other read()/edit() paths still need their normal semantic regression checks when core is changed.
