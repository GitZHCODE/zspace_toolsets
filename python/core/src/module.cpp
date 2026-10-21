#include <pybind11/pybind11.h>

#include "bindings.h"
#include <zspace/zInterface/objects/zObjMeshField.h>

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
  m.doc() = "zSpace Core SDK bindings: zVector, zColor, zObjMesh/zFnMesh, "
            "zObjGraph/zFnGraph.";

  zspace_python::bind_core_types(m);
  zspace_python::bind_obj_mesh(m);
  zspace_python::bind_fn_mesh(m);
  zspace_python::bind_obj_graph(m);
  zspace_python::bind_fn_graph(m);
  py::class_<zSpace::zObjMeshScalarField>(m, "ObjMeshScalarField").def(py::init<>());
}
