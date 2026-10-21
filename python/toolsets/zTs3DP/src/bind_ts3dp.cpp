// Bindings for zTsSlicing (zTsSlicing.h) - one slicer operator (TopologySlicer topology slicing,
// or HarmonicSlicer harmonic slicing - see setSlicerMode()); each operator
// post-processes its own result into print paths with the zPrintSyn functions
// (toolpath extraction, cleanup, projection, chaining).
//
// computeSlicer() then computePostProcess() produces getResultantPrintPath() - a
// single continuous bottom-to-top PrintPath chained across every layer. That
// is normally the one result you want out of this whole pipeline.
//
// Only one operator exists: the default (Topology) one from the constructor,
// replaced by a new one on every setSlicerMode() call (old results are
// discarded). getSettings() (one SlicerSettings for both) is held by zTsSlicing and
// editable at any time; getWidthGroupCount()/getIsoMeshes()/... read
// the Topology operator and are 1/empty when it is not the one instantiated.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <zspace/zToolsets/fabrication/slicing/zTsSlicing.h>

namespace py = pybind11;
using namespace zSpace;

namespace zspace_python {

void bind_ts3dp(py::module_ &m) {
  py::class_<zTsSlicing>(m, "zTsSlicing")
      .def(py::init<>())
      .def_static("supportsSlicerMode", &zTsSlicing::supportsSlicerMode)

      .def("setInputMesh", &zTsSlicing::setInputMesh, py::arg("mesh"),
           py::arg("face_tags") = zIntArray(),
           py::arg("face_width_ids") = zIntArray())
      .def("setBracingGraph", &zTsSlicing::setBracingGraph, py::arg("graph"))
      .def("setSlotAnchor", &zTsSlicing::setSlotAnchor, py::arg("anchor"))
      .def("clearSlotAnchor", &zTsSlicing::clearSlotAnchor)
      .def("hasSlotAnchor", &zTsSlicing::hasSlotAnchor)
      .def("getSlotAnchor", &zTsSlicing::getSlotAnchor)
      .def("getFaceTags", &zTsSlicing::getFaceTags)
      .def("getFaceWidthIds", &zTsSlicing::getFaceWidthIds)
      .def("setSlicerMode", &zTsSlicing::setSlicerMode, py::arg("mode"))
      .def("getSlicerMode", &zTsSlicing::getSlicerMode)

      .def("computeSlicer", &zTsSlicing::computeSlicer, py::arg("use_parallel") = true, py::call_guard<py::gil_scoped_release>())
      .def("computePrintPath", &zTsSlicing::computePrintPath, py::arg("use_parallel") = true, py::call_guard<py::gil_scoped_release>())
      .def("computeBasicStructure", &zTsSlicing::computeBasicStructure, py::arg("use_parallel") = true, py::call_guard<py::gil_scoped_release>())
      .def("computePostProcess", &zTsSlicing::computePostProcess, py::arg("use_parallel") = true, py::call_guard<py::gil_scoped_release>())

      .def("isLoaded", &zTsSlicing::isLoaded)
      .def("getStatus", &zTsSlicing::getStatus)

      .def("getSettings",
           static_cast<SlicerSettings &(zTsSlicing::*)()>(&zTsSlicing::getSettings),
           py::return_value_policy::reference_internal)
      .def("getWidthGroupCount", &zTsSlicing::getWidthGroupCount)
      .def("getPrintSynSetting",
           static_cast<PrintSynSetting &(zTsSlicing::*)()>(
               &zTsSlicing::getPrintSynSetting),
           py::return_value_policy::reference_internal)
      // Read-only views of the selected operator's computePostProcess()
      // results (copied to Python lists / bound to this zTsSlicing).
      .def("get2DPrintPaths", &zTsSlicing::get2DPrintPaths)
      .def("get3DPrintPaths", &zTsSlicing::get3DPrintPaths)
      .def("getPrintHeightDomain", [](const zTsSlicing &s) {
        const auto &d = s.getPrintHeightDomain(); return py::make_tuple(d.min, d.max);
      })
      .def("getPrintWidthDomain", [](const zTsSlicing &s) {
        const auto &d = s.getPrintWidthDomain(); return py::make_tuple(d.min, d.max);
      })
      .def("getPrintLength", &zTsSlicing::getPrintLength)
      .def("getPrintVolume", &zTsSlicing::getPrintVolume)
      .def("hasPrintStatistics", &zTsSlicing::hasPrintStatistics)
      .def("getResultantPrintPath", &zTsSlicing::getResultantPrintPath,
           py::return_value_policy::reference_internal)
      .def("getIsoMeshes",
           static_cast<zObjMeshArray &(zTsSlicing::*)()>(&zTsSlicing::getIsoMeshes))
      .def("getIsoMeshes2D",
           static_cast<zObjMeshArray &(zTsSlicing::*)()>(
               &zTsSlicing::getIsoMeshes2D))
      .def("getScalarField",
           static_cast<zObjMeshScalarField &(zTsSlicing::*)()>(
               &zTsSlicing::getScalarField),
           py::return_value_policy::reference_internal)
      .def("getPathFields",
           static_cast<std::vector<zScalarArray> &(zTsSlicing::*)()>(
               &zTsSlicing::getPathFields))
      .def("getWidthFields",
           static_cast<std::vector<zScalarArray> &(zTsSlicing::*)()>(
               &zTsSlicing::getWidthFields))
      .def("getTrackFields",
           static_cast<std::vector<zScalarArray> &(zTsSlicing::*)()>(
               &zTsSlicing::getTrackFields))
      .def("getInputMesh",
           static_cast<zObjMesh &(zTsSlicing::*)()>(&zTsSlicing::getInputMesh),
           py::return_value_policy::reference_internal)
      .def("getBracingGraph",
           static_cast<zObjGraph &(zTsSlicing::*)()>(&zTsSlicing::getBracingGraph),
           py::return_value_policy::reference_internal)
      ;
}

}  // namespace zspace_python
