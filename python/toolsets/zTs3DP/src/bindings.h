#pragma once

#include <pybind11/pybind11.h>

namespace zspace_python {

// Order matters: bind_types (SlicerMode, SlicerSettings,
// PrintSynSetting, ParameterGraph, PrintPath,
// ObjMeshScalarField) must run before bind_ts3dp, which uses them as
// parameters/return values.
//
// TopologySlicer, zPrintSyn and zTsUtils are intentionally not bound to Python -
// zTs3DP is the only exposed entry point into the slicer toolset.
void bind_types(pybind11::module_ &m);
void bind_ts3dp(pybind11::module_ &m);      // zTs3DP

}  // namespace zspace_python
