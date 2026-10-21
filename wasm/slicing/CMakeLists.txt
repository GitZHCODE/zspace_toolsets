cmake_minimum_required(VERSION 3.21)
project(zspace_slicing_wasm LANGUAGES CXX)
if(NOT EMSCRIPTEN)
  message(FATAL_ERROR "Configure with Emscripten/emcmake")
endif()
get_filename_component(TOOLSETS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(ZSPACE_CORE_DIR "${TOOLSETS_ROOT}/../zspace_core" CACHE PATH "Core source checkout")
set(ZSPACE_TETGEN_SOURCE_DIR "${TOOLSETS_ROOT}/wasm/cache/tetgen" CACHE PATH "TetGen source checkout")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS
  "${ZSPACE_CORE_DIR}/src/zCore/*.cpp"
  "${ZSPACE_CORE_DIR}/src/zInterface/*.cpp"
  "${ZSPACE_CORE_DIR}/src/zIO/*.cpp")
list(FILTER sources EXCLUDE REGEX "/(zFnComputeMesh|zObjectComputeMesh|zObjectComputeField2D)\\.cpp$")
list(FILTER sources EXCLUDE REGEX "/(zFnMeshField|zFnPointField)\\.cpp$")
list(APPEND sources "${CMAKE_CURRENT_SOURCE_DIR}/mesh_field_compat.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/point_field_compat.cpp")
set(harmonic OFF)
if(EXISTS "${ZSPACE_TETGEN_SOURCE_DIR}/tetgen.cxx" AND EXISTS "${ZSPACE_TETGEN_SOURCE_DIR}/predicates.cxx")
  set(harmonic ON)
  list(APPEND sources "${ZSPACE_TETGEN_SOURCE_DIR}/tetgen.cxx" "${ZSPACE_TETGEN_SOURCE_DIR}/predicates.cxx")
endif()
include("${TOOLSETS_ROOT}/cmake/SlicingSources.cmake")
zspace_slicing_sources(slicing_sources "${TOOLSETS_ROOT}" "${harmonic}")
# Browser slicing is compute-only; USD/file IO belongs to core.
list(FILTER slicing_sources EXCLUDE REGEX "/utils/IO\\.cpp$")
list(APPEND sources ${slicing_sources})
add_executable(zspace_slicing ${sources} "${CMAKE_CURRENT_SOURCE_DIR}/zspace_slicing_bridge.cpp")
set_target_properties(zspace_slicing PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES
  SUFFIX ".js" RUNTIME_OUTPUT_DIRECTORY "${TOOLSETS_ROOT}/wasm/out")
target_include_directories(zspace_slicing PRIVATE "${ZSPACE_CORE_DIR}/wasm/bridge" "${TOOLSETS_ROOT}/include"
  "${TOOLSETS_ROOT}/wasm/bridge" "${ZSPACE_CORE_DIR}" "${ZSPACE_CORE_DIR}/include"
  "${ZSPACE_CORE_DIR}/third_party" "${ZSPACE_CORE_DIR}/third_party/depends")
target_compile_definitions(zspace_slicing PRIVATE ZSPACE_STATIC_LIBRARY ZSPACE_TOOLSETS_STATIC_LIBRARY EIGEN_DONT_PARALLELIZE)
if(harmonic)
  target_include_directories(zspace_slicing PRIVATE "${ZSPACE_TETGEN_SOURCE_DIR}")
  target_compile_definitions(zspace_slicing PRIVATE ZSPACE_TOOLSETS_WITH_TETGEN TETLIBRARY)
endif()
target_compile_options(zspace_slicing PRIVATE -O2 -fexceptions -Uunix -include "${ZSPACE_CORE_DIR}/wasm/bridge/zspace_wasm_compat.h")
target_precompile_headers(zspace_slicing PRIVATE <depends/Eigen/Dense> <nlohmann/json.hpp>)
target_link_options(zspace_slicing PRIVATE -O2 -fexceptions "-sMODULARIZE=1" "-sEXPORT_ES6=1"
  "-sENVIRONMENT=web,worker,node" "-sALLOW_MEMORY_GROWTH=1" "-sSTACK_SIZE=2097152"
  "-sEXPORTED_RUNTIME_METHODS=['HEAPU32','UTF8ToString','ccall']")
message(STATUS "Slicing WASM harmonic support: ${harmonic}")
