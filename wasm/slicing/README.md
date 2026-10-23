# Combined Toolsets WASM runtime

Build with `../scripts/build-slicing-wasm.ps1`. Outputs are
`../out/zspace_toolset.js` and `../out/zspace_toolset.wasm`.

Core is included with add_subdirectory and linked as reusable objects. Core owns
its source list, geometry bindings, field compatibility shims, USD SDK resolution
and compilation settings. Toolsets owns its algorithms and slicer bridge only.
Allocation and Core geometry C ABI exports come from Core in the same memory.

Slicer array ABI 2 replaces the old binary packet entrypoints. Inputs and results
use numeric array pointers/counts, without JSON intermediaries. Core pointer/count
exports remain available in the same runtime.
See [array ABI](../../docs/architecture/slicing-wasm-arrays.md).
The unrelated zspace_toolsets runtime (plural) remains separate.

Thread pool and memory settings belong to the final runtime. Core resolves SDK
paths. No additional standalone Core WASM is needed for this runtime's Core APIs.

The shared WASM compilation and optimized linking use -O2, without explicit
NDEBUG, SIMD or LTO flags. Keep ZSPACE_WASM_FAST_LINK disabled for optimized
builds. Imported prebuilt SDK archives retain their original compilation options.
