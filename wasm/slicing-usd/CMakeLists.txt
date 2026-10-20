cmake_minimum_required(VERSION 3.24)
project(zspace_slicing_usd LANGUAGES CXX)
get_filename_component(TOOLSETS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(ZSPACE_CORE_DIR "${TOOLSETS_ROOT}/../zspace_core" CACHE PATH "Core source checkout")
add_subdirectory("${ZSPACE_CORE_DIR}/wasm/openusd" core)
target_sources(zspace_core PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/slicing_usd_bridge.cpp"
  "${TOOLSETS_ROOT}/src/zToolsets/fabrication/slicing/utils/IO.cpp"
  "${TOOLSETS_ROOT}/src/zToolsets/fabrication/slicing/utils/Utils.cpp")
target_include_directories(zspace_core PRIVATE "${TOOLSETS_ROOT}/include" "${TOOLSETS_ROOT}/wasm/bridge")
target_compile_definitions(zspace_core PRIVATE ZSPACE_TOOLSETS_STATIC_LIBRARY ZSPACE_TS3DP_OPENUSD)
target_compile_options(zspace_core PRIVATE -Uunix)
set_target_properties(zspace_core PROPERTIES OUTPUT_NAME zspace_slicing_usd
  RUNTIME_OUTPUT_DIRECTORY "${TOOLSETS_ROOT}/wasm/out/slicing-usd")

# Flatten the SDK archives exactly once. Transitive CMake static-library ordering
# repeats usd/usdGeom, which produces duplicate definitions inside whole-archive.
# Whole archives preserve OpenUSD's static schema/plugin registrations.
file(GLOB usdArchives "${ZSPACE_OPENUSD_ROOT}/lib/libusd_*.a")
if(NOT usdArchives)
  message(FATAL_ERROR "No static OpenUSD archives in ${ZSPACE_OPENUSD_ROOT}/lib")
endif()
set_property(TARGET zspace_core PROPERTY LINK_LIBRARIES "")
foreach(archive IN LISTS usdArchives)
  target_link_libraries(zspace_core PRIVATE "$<LINK_LIBRARY:WHOLE_ARCHIVE,${archive}>")
endforeach()
target_link_libraries(zspace_core PRIVATE "${ZSPACE_OPENUSD_ROOT}/lib/libtbb.a")
target_include_directories(zspace_core SYSTEM PRIVATE "${ZSPACE_OPENUSD_ROOT}/include")
target_compile_definitions(zspace_core PRIVATE PXR_STATIC)
# The archive target's transitive embed flags were removed with its dependency
# graph above. Ship the schema/plugin resources at the SDK's runtime /usd path.
target_link_options(zspace_core PRIVATE
  "SHELL:--embed-file \"${ZSPACE_OPENUSD_ROOT}/lib/usd@/usd\"")