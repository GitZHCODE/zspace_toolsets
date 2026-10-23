#include <pybind11/pybind11.h>

#include "bindings.h"
#include <zspace/zToolsets/fabrication/slicing/zTsSlicing.h>
#if defined(ZSPACE_TOOLSETS_WITH_USD_IO)
#include <zspace/zToolsets/fabrication/slicing/utils/IO.h>
#endif
#include <pybind11/stl.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

namespace py = pybind11;

PYBIND11_MODULE(_slicer, m) {
  py::module_::import("zspace.core");
  m.doc() = "zspace.toolset.slicer: pure geometry pipeline and settings. TopologySlicer, "
            "zPrintSyn and zTsUtils are intentionally not exposed to Python.";

  zspace_python::bind_types(m);
  zspace_python::bind_ts3dp(m);
  m.def("buildPrintMesh", [](const zSpace::zTsSlicing &solver) {
    auto mesh = std::make_unique<zSpace::zObjMesh>();
    std::string error;
    if (!zSpace::zTsUtils::buildPrintMesh(solver.get3DPrintPaths(), *mesh, &error))
      throw std::runtime_error(error);
    return mesh;
  });
#if defined(ZSPACE_TOOLSETS_WITH_USD_IO)
  using namespace zSpace;
  py::class_<zTsIO::SlicerUSDInput>(m, "SlicerUSDInput")
      .def_readonly("Block", &zTsIO::SlicerUSDInput::Block)
      .def_readonly("Bracing", &zTsIO::SlicerUSDInput::Bracing)
      .def_readonly("FaceTags", &zTsIO::SlicerUSDInput::FaceTags)
      .def_readonly("FaceWidthIds", &zTsIO::SlicerUSDInput::FaceWidthIds)
      .def_readonly("StartPoint", &zTsIO::SlicerUSDInput::StartPoint)
      .def_readonly("Settings", &zTsIO::SlicerUSDInput::Settings)
      .def_readonly("PrintSettings", &zTsIO::SlicerUSDInput::PrintSettings)
      .def_readonly("BracingPrintWidths", &zTsIO::SlicerUSDInput::BracingPrintWidths)
      .def_readonly("UnknownSettings", &zTsIO::SlicerUSDInput::UnknownSettings)
      .def_readonly("MetersPerUnit", &zTsIO::SlicerUSDInput::MetersPerUnit)
      .def_readonly("UpAxis", &zTsIO::SlicerUSDInput::UpAxis);
  m.def("importUSD", [](const std::string &path) {
    auto input = std::make_unique<zTsIO::SlicerUSDInput>(); std::string error;
    if (!zTsIO::importUSD(path, *input, &error)) throw std::runtime_error(error);
    return input;
  });
  m.def("exportUSD", [](const PrintPath &path, const std::string &input, const std::string &output, int layers, int blockid, const zTsSlicing *solver) {
    std::string error;
    if (!zTsIO::exportUSD(path,input,output,layers,&error,blockid,solver ? &solver->get3DPrintPaths() : nullptr)) throw std::runtime_error(error);
  }, py::arg("path"),py::arg("input"),py::arg("output"),py::arg("layers"),py::arg("blockid")=0,py::arg("solver")=nullptr);
#endif
  m.attr("zTs3DP") = m.attr("zTsSlicing");
}
