# Native geometry bindings linked against the headless core source build.
pybind11_add_module(_core
    src/module.cpp
    src/bind_core_types.cpp
    src/bind_obj_mesh.cpp
    src/bind_fn_mesh.cpp
    src/bind_obj_graph.cpp
    src/bind_fn_graph.cpp
)

target_link_libraries(_core PRIVATE zspace::core zspace::interface)
set(ZSPACE_PYTHON_PACKAGE_DIR "${ZSPACE_PYTHON_PACKAGE_DIR}/core")
# Plain LIBRARY_/RUNTIME_OUTPUT_DIRECTORY get a per-config subdirectory
# appended by multi-config generators (VS: Release/, Debug/, ...) unless the
# per-config properties are set explicitly, so set both.
set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${ZSPACE_PYTHON_PACKAGE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${ZSPACE_PYTHON_PACKAGE_DIR}"
)
foreach(_cfg IN ITEMS DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(_core PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY_${_cfg} "${ZSPACE_PYTHON_PACKAGE_DIR}"
        RUNTIME_OUTPUT_DIRECTORY_${_cfg} "${ZSPACE_PYTHON_PACKAGE_DIR}"
    )
endforeach()

# The zSpace Core DLLs (zSpace_Core.dll, zSpace_Interface.dll, ...) must sit
# next to _core.pyd for Python to load it at import time.
add_custom_command(TARGET _core POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:zSpace_Core>" "$<TARGET_FILE:zSpace_Interface>"
    "${ZSPACE_PYTHON_PACKAGE_DIR}" VERBATIM)
