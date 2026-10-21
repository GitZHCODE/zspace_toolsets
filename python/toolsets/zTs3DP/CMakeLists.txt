# Bind the v3 library; do not include implementation .cpp files in the bindings.
pybind11_add_module(_slicer src/module.cpp src/bind_types.cpp src/bind_ts3dp.cpp)
target_link_libraries(_slicer PRIVATE zspace::toolsets)
if(ZSPACE_PYTHON_WITH_USD)
  target_compile_definitions(_slicer PRIVATE ZSPACE_TS3DP_OPENUSD)
endif()
set(_slicer_package "${ZSPACE_PYTHON_PACKAGE_DIR}/toolset/slicer")
set_target_properties(_slicer PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${_slicer_package}"
  RUNTIME_OUTPUT_DIRECTORY "${_slicer_package}")
foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
  set_target_properties(_slicer PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY_${_cfg} "${_slicer_package}"
    RUNTIME_OUTPUT_DIRECTORY_${_cfg} "${_slicer_package}")
endforeach()
add_custom_command(TARGET _slicer POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "$<TARGET_FILE:zSpace_Toolsets>" "$<TARGET_FILE:zSpace_IO>"
    "${_slicer_package}" VERBATIM)
