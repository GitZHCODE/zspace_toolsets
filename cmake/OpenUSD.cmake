# Source-built OpenUSD 25.08, Release x64, monolithic, without Python/imaging.
if(NOT TARGET zspace_toolsets_usd)
  set(ZSPACE_TOOLSETS_OPENUSD_ROOT "${CMAKE_CURRENT_LIST_DIR}/../cpp/depends/openusd-minimal" CACHE PATH "Source-built minimal OpenUSD SDK")
  # Migrate the removed downloaded SDK path from existing alice2 CMake caches.
  if(ZSPACE_TOOLSETS_OPENUSD_ROOT MATCHES "openusd/25[.]08/?$")
    set(ZSPACE_TOOLSETS_OPENUSD_ROOT "${CMAKE_CURRENT_LIST_DIR}/../cpp/depends/openusd-minimal" CACHE PATH "Source-built minimal OpenUSD SDK" FORCE)
    unset(pxr_DIR CACHE)
  endif()
  get_filename_component(_usd_root "${ZSPACE_TOOLSETS_OPENUSD_ROOT}" ABSOLUTE)
  if(NOT EXISTS "${_usd_root}/pxrConfig.cmake")
    message(FATAL_ERROR "Build the minimal SDK first: missing ${_usd_root}/pxrConfig.cmake")
  endif()
  find_package(pxr CONFIG REQUIRED PATHS "${_usd_root}" NO_DEFAULT_PATH)
  if(NOT TARGET usd_ms)
    message(FATAL_ERROR "The toolsets USD SDK must be built with --build-monolithic")
  endif()
  if(NOT TARGET TBB::tbb)
    message(FATAL_ERROR "The OpenUSD SDK did not provide its TBB dependency")
  endif()
  if(NOT EXISTS "${_usd_root}/bin/tbb.dll")
    message(FATAL_ERROR "Missing OpenUSD runtime dependency: ${_usd_root}/bin/tbb.dll")
  endif()
  # This SDK and its TBB were built in Release, including for RelWithDebInfo consumers.
  set_property(TARGET TBB::tbb PROPERTY INTERFACE_COMPILE_DEFINITIONS "")
  add_library(zspace_toolsets_usd INTERFACE)
  set_property(TARGET zspace_toolsets_usd PROPERTY ZSPACE_USD_RESOURCE_DIR "${_usd_root}/lib/usd")
  target_compile_definitions(zspace_toolsets_usd INTERFACE ZSPACE_TS3DP_OPENUSD)
  target_link_libraries(zspace_toolsets_usd INTERFACE usd_ms)
  target_include_directories(zspace_toolsets_usd BEFORE INTERFACE "${_usd_root}/include")
endif()

# Call for an executable or pybind MODULE/SHARED target that consumes USD.
# Plain link signature matches alice2's existing target_link_libraries calls.
function(zspace_toolsets_enable_usd target)
  if(NOT TARGET "${target}")
    message(FATAL_ERROR "Unknown USD consumer target: ${target}")
  endif()
  get_target_property(_usd_attached "${target}" ZSPACE_TOOLSETS_USD_ATTACHED)
  if(_usd_attached)
    return()
  endif()
  target_link_libraries("${target}" zspace_toolsets_usd)
  get_target_property(_usd_resources zspace_toolsets_usd ZSPACE_USD_RESOURCE_DIR)
  if(WIN32)
    add_custom_command(TARGET "${target}" POST_BUILD
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "$<TARGET_FILE:usd_ms>" "$<TARGET_FILE_DIR:${target}>"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "$<TARGET_FILE:TBB::tbb>" "$<TARGET_FILE_DIR:${target}>"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
        "${_usd_resources}" "$<TARGET_FILE_DIR:${target}>/usd"
      COMMENT "Deploying OpenUSD, TBB and USD schema resources"
      VERBATIM)
  endif()
  set_property(TARGET "${target}" PROPERTY ZSPACE_TOOLSETS_USD_ATTACHED TRUE)
endfunction()
