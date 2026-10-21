# Slicer implementation

The earlier 3DP prototype and its pseudocode have been retired.
The maintained implementation is `fabrication/slicing/zTsSlicing`.

See [tagged slicing](tagged-slicing.md) and
[the current example](../../examples/3dp/example_3dp_sdf_slicing.cpp).

Topology and Harmonic share field utilities and PrintSyn. Call `computeSlicer()`
followed by `computePostProcess()` for the final connected path, or
`computePrintPath()` for per-layer extracted and mapped paths without PrintSyn.
