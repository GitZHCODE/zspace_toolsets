#include <pybind11/pybind11.h>

#include "bindings.h"
#include <pybind11/stl.h>
#include <zspace/zToolsets/fabrication/slicing/utils/Utils.h>

namespace py = pybind11;

PYBIND11_MODULE(_slicer, m) {
  py::module_::import("zspace.core");
  m.doc() = "zspace.toolset.slicer: pipeline and settings with optional USD IO. TopologySlicer, "
            "zPrintSyn and zTsUtils are intentionally not exposed to Python.";

  zspace_python::bind_types(m);
  zspace_python::bind_ts3dp(m);
  m.attr("zTs3DP") = m.attr("zTsSlicing");
#if defined(ZSPACE_TS3DP_OPENUSD)
  m.attr("HAS_USD") = true;
  using zSpace::zTsIO::SlicerUSDInput;
  py::class_<SlicerUSDInput>(m, "SlicerUSDInput")
      .def_readonly("block", &SlicerUSDInput::Block)
      .def_readonly("face_tags", &SlicerUSDInput::FaceTags)
      .def_readonly("face_width_ids", &SlicerUSDInput::FaceWidthIds)
      .def_readonly("bracing_curves", &SlicerUSDInput::BracingCurves)
      .def_readonly("bracing_width_ids", &SlicerUSDInput::BracingWidthIds)
      .def_readonly("bracing_print_widths", &SlicerUSDInput::BracingPrintWidths)
      .def_readonly("slotpoint", &SlicerUSDInput::StartPoint)
      .def_readonly("settings_json", &SlicerUSDInput::SettingsJSON)
      .def_readonly("meters_per_unit", &SlicerUSDInput::MetersPerUnit)
      .def_readonly("up_axis", &SlicerUSDInput::UpAxis);
  m.def("read_input_usd", [](const std::string &file) {
    SlicerUSDInput input;
    std::string error;
    if (!zSpace::zTsIO::readSlicerInputUSD(file, input, &error))
      throw std::runtime_error(error);
    return input;
  }, py::arg("file"));
  m.def("export_printpath_usd", [](const zSpace::PrintPath &path,
      const std::string &input, const std::string &output, int layerCount, int blockid) {
    std::string error;
    if (!zSpace::zTsIO::exportPrintPathUSD(path, input, output, layerCount, &error, blockid))
      throw std::runtime_error(error);
  }, py::arg("path"), py::arg("input_file"), py::arg("output_file"), py::arg("layercount"), py::arg("blockid") = 0);
#else
  m.attr("HAS_USD") = false;
#endif
}
