# Toolsets browser runtime

Build with `scripts/build-wasm.bat` from this folder, or
`wasm/scripts/build-wasm.bat` from the repository root. `scripts/build-wasm.ps1`
accepts `-CoreDir` and `-EmsdkDir`. Generated build/cache/tmp/out files stay ignored.

The toolsets bridge compiles against the current sibling core sources. Its API
and heap are independent of the core viewer runtime. Follow
[typed transfer contract](../docs/architecture/typed-wasm-transfers.md).

Preserve existing C exports; add transfer APIs without changing the legacy ABI.
Refresh the viewer's packaged JS/WASM pair together only when runtime integration
is in scope. Read applicable viewer instructions before editing adapters.

## Slicing OpenMP runtime

Build `wasm/scripts/build-slicing-wasm.ps1 -EmsdkDir <emsdk>` to produce
`wasm/out/zspace_slicing.js` and `zspace_slicing.wasm`. Both files must be deployed
together. This compute-only target compiles sibling core and slicer sources;
it does not link OpenUSD.

The target enables `-pthread -fopenmp` at compile and link time. It pre-creates
`navigator.hardwareConcurrency` minus one pthread workers (the calling Worker is also a
compute thread), limits OpenMP to that actual pool capacity, disables dynamic/nested
teams, and keeps Eigen's internal parallelism disabled. The transfer input
`useParallel` defaults to `false`; `true` enables the existing slicer and PrintSyn
OpenMP loops. Results report `useParallel` and the observed `threadCount`.
Serial and parallel execution share the same threaded binary.

The viewer runs this module inside its slicing Worker, leaving the UI responsive.
Its **Compute** action selects serial execution; **Parallel compute** selects
OpenMP execution. One request runs at a time per session. Cancelling or disposing
the session terminates the owning Worker and its child workers. Node consumers
must call `module.PThread.terminateAllThreads()` when finished.

Both modes require SharedArrayBuffer and a secure, cross-origin-isolated browser
context (localhost or HTTPS, COOP `same-origin`, COEP `require-corp`). The viewer's
Vite development and preview servers already provide these headers. Static
production hosting must provide them too. The OpenUSD runtime has a separate
thread pool; it is not part of the slicer's OpenMP team.

Thread capacity is selected at module initialization from the browser-reported
logical processor count, falling back to one when unavailable. It is not live
CPU-load detection; browser privacy limits may lower the reported count. Recreate
the session/reload the page to pick up a changed capacity. More threads can cost
more memory and do not guarantee a speedup. OpenMP dynamic team adjustment stays
disabled so the parallel action requests the full available team.
