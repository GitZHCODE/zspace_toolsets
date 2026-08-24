# zspace_toolsets

Modern zSpace toolsets are built as focused libraries on top of `zspace_core`.

The pre-v3 header/source/project layout has been moved to `legacy/pre-v3`. New work should use:

```text
include/zspace/zToolsets/
src/zToolsets/
examples/
tests/
```

The first migrated toolset is `zTs3DP`, a 3D-printing pipeline for SDF-based slicing, unrolled slice fields, contour extraction, and print-path synthesis.

Build:

```bat
scripts\build_toolsets.bat
```
