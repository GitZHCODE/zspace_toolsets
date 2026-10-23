# Bind the v3 library; do not include implementation .cpp files in the bindings.
pybind11_add_module(_slicer src/module.cpp src/bind_types.cpp src/bind_ts3dp.cpp)
target_link_libraries(_slicer PRIVATE zspace::toolsets)
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
    "$<TARGET_FILE:zSpace_Toolsets>"
    "${ZSPACE_PYTHON_PACKAGE_DIR}/runtime" VERBATIM)

target_include_directories(_slicer PRIVATE "${ZSPACE_CORE_SOURCE}/python/core/src")

# Core's binding target deploys IO and USD once into the shared runtime.
add_dependencies(_slicer _core)

if(WIN32 AND ZSPACE_TOOLSETS_SLICING_OPENMP)
  set(ZSPACE_PYTHON_OPENMP_RUNTIME "" CACHE FILEPATH "MSVC OpenMP runtime DLL")
  if(NOT ZSPACE_PYTHON_OPENMP_RUNTIME)
    get_filename_component(_vc_root "${CMAKE_CXX_COMPILER}/../../../../../../.." ABSOLUTE)
    file(GLOB _openmp_candidates "${_vc_root}/Redist/MSVC/*/x64/Microsoft.VC*.OpenMP/vcomp140.dll")
    list(SORT _openmp_candidates COMPARE NATURAL ORDER DESCENDING)
    if(_openmp_candidates)
      list(GET _openmp_candidates 0 ZSPACE_PYTHON_OPENMP_RUNTIME)
    endif()
  endif()
  if(NOT EXISTS "${ZSPACE_PYTHON_OPENMP_RUNTIME}")
    message(FATAL_ERROR "Set ZSPACE_PYTHON_OPENMP_RUNTIME to vcomp140.dll")
  endif()
  add_custom_command(TARGET _slicer POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy_if_different "${ZSPACE_PYTHON_OPENMP_RUNTIME}" "${ZSPACE_PYTHON_PACKAGE_DIR}/runtime" VERBATIM)
endif()

