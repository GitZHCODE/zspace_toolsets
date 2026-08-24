cmake_minimum_required(VERSION 3.21)

project(zspace_toolsets VERSION 3.0.0 LANGUAGES CXX)

option(ZSPACE_TOOLSETS_BUILD_SHARED "Build zSpace toolsets as a DLL" ON)
option(ZSPACE_TOOLSETS_BUILD_TESTS "Build zSpace toolsets smoke tests" OFF)
option(ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE "Build against a sibling zspace_core checkout" ON)
set(ZSPACE_CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../zspace_core" CACHE PATH "Path to local zspace_core checkout")
set(ZSPACE_SDK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/sdk/zspace" CACHE PATH "Path to installed zSpace SDK")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(ZSPACE_TOOLSETS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

if(ZSPACE_TOOLSETS_BUILD_SHARED)
  set(ZSPACE_TOOLSETS_LIBRARY_TYPE SHARED)
else()
  set(ZSPACE_TOOLSETS_LIBRARY_TYPE STATIC)
endif()

if(TARGET zspace::core AND TARGET zspace::interface AND TARGET zspace::io)
  message(STATUS "zspace_toolsets: using existing zSpace targets from parent project")
elseif(ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE)
  if(NOT EXISTS "${ZSPACE_CORE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE=ON requires ZSPACE_CORE_DIR to point to zspace_core")
  endif()

  set(ZSPACE_BUILD_INTERFACE ON CACHE BOOL "Build zSpace interface DLL" FORCE)
  set(ZSPACE_BUILD_IO ON CACHE BOOL "Build zSpace IO DLL" FORCE)
  set(ZSPACE_BUILD_DISPLAY OFF CACHE BOOL "Build zSpace display DLL" FORCE)
  set(ZSPACE_BUILD_INTEROP OFF CACHE BOOL "Build zSpace interop DLL" FORCE)
  set(ZSPACE_BUILD_TESTS OFF CACHE BOOL "Build zSpace tests" FORCE)
  set(ZSPACE_BUILD_SHARED ON CACHE BOOL "Build zSpace DLLs" FORCE)
  set(ZSPACE_WITH_OPENGL OFF CACHE BOOL "Disable zSpace display OpenGL" FORCE)
  set(ZSPACE_IO_WITH_TINYUSDZ OFF CACHE BOOL "Disable TinyUSDZ download for toolsets scaffold" FORCE)

  add_subdirectory("${ZSPACE_CORE_DIR}" "${CMAKE_BINARY_DIR}/zspace_core_build")
else()
  list(PREPEND CMAKE_PREFIX_PATH "${ZSPACE_SDK_DIR}")
  find_package(zspace 3 CONFIG REQUIRED)
endif()

file(GLOB_RECURSE ZSPACE_TOOLSETS_SOURCES CONFIGURE_DEPENDS
  "${ZSPACE_TOOLSETS_ROOT}/src/zToolsets/*.cpp"
)
file(GLOB_RECURSE ZSPACE_TOOLSETS_HEADERS CONFIGURE_DEPENDS
  "${ZSPACE_TOOLSETS_ROOT}/include/zspace/*.h"
)

add_library(zSpace_Toolsets ${ZSPACE_TOOLSETS_LIBRARY_TYPE}
  ${ZSPACE_TOOLSETS_SOURCES}
  ${ZSPACE_TOOLSETS_HEADERS}
)
add_library(zspace::toolsets ALIAS zSpace_Toolsets)
set_target_properties(zSpace_Toolsets PROPERTIES EXPORT_NAME toolsets)

target_include_directories(zSpace_Toolsets
  PUBLIC
    $<BUILD_INTERFACE:${ZSPACE_TOOLSETS_ROOT}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(zSpace_Toolsets PUBLIC zspace::core zspace::interface zspace::io)

target_compile_definitions(zSpace_Toolsets
  PUBLIC
    $<$<BOOL:${ZSPACE_TOOLSETS_BUILD_SHARED}>:ZSPACE_TOOLSETS_DYNAMIC_LIBRARY>
    $<$<NOT:$<BOOL:${ZSPACE_TOOLSETS_BUILD_SHARED}>>:ZSPACE_TOOLSETS_STATIC_LIBRARY>
  PRIVATE
    $<$<BOOL:${ZSPACE_TOOLSETS_BUILD_SHARED}>:ZSPACE_TOOLSETS_EXPORTS>
)

if(MSVC)
  target_compile_options(zSpace_Toolsets PRIVATE /MP /W3 /permissive- /wd4251)
  target_compile_definitions(zSpace_Toolsets PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

install(TARGETS zSpace_Toolsets
  EXPORT zspaceToolsetsTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION bin
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include
)
install(DIRECTORY "${ZSPACE_TOOLSETS_ROOT}/include/zspace" DESTINATION include)
install(EXPORT zspaceToolsetsTargets FILE zspaceToolsetsTargets.cmake NAMESPACE zspace:: DESTINATION lib/cmake/zspace_toolsets)

if(ZSPACE_TOOLSETS_BUILD_TESTS)
  include(CTest)
  add_executable(zspace_toolsets_3dp_smoke tests/smoke/zspace_toolsets_3dp_smoke.cpp)
  target_link_libraries(zspace_toolsets_3dp_smoke PRIVATE zspace::toolsets)
  add_test(NAME zspace_toolsets.3dp.smoke COMMAND zspace_toolsets_3dp_smoke)
endif()

message(STATUS "zspace_toolsets Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  zSpace source mode: ${ZSPACE_TOOLSETS_USE_ZSPACE_SOURCE}")
message(STATUS "  zSpace core dir: ${ZSPACE_CORE_DIR}")

