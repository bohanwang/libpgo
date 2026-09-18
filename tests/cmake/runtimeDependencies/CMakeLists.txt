cmake_minimum_required(VERSION 3.28 FATAL_ERROR)

project(runtimeDependenciesSmoke LANGUAGES CXX)

if(NOT LIBPGO_SOURCE_DIR)
  message(FATAL_ERROR "LIBPGO_SOURCE_DIR is required")
endif()

include("${LIBPGO_SOURCE_DIR}/CMakeModules/runtimeDependencies.cmake")
pgo_configure_runtime_dependencies()

add_library(runtime_probe MODULE probe.cpp)
pgo_attach_runtime_dependencies(runtime_probe WHEEL_EXTENSION)

if(APPLE)
  get_target_property(_pgo_probe_rpath runtime_probe BUILD_RPATH)
  if(NOT _pgo_probe_rpath STREQUAL "@loader_path/../..")
    message(FATAL_ERROR "Unexpected macOS wheel RPATH: ${_pgo_probe_rpath}")
  endif()
elseif(UNIX)
  get_target_property(_pgo_probe_rpath runtime_probe BUILD_RPATH)
  if(NOT _pgo_probe_rpath STREQUAL "$ORIGIN/../..")
    message(FATAL_ERROR "Unexpected Linux wheel RPATH: ${_pgo_probe_rpath}")
  endif()
elseif(WIN32 AND TARGET pgo_windows_runtime)
  message(FATAL_ERROR "WHEEL layout must not create source DLL staging")
endif()
