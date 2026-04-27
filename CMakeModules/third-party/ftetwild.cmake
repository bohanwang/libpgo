message(STATUS "Loading fTetWild...")

include(FetchContent)

set(_PGO_SAVED_CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD}")
set(_PGO_SAVED_CMAKE_CXX_STANDARD_REQUIRED "${CMAKE_CXX_STANDARD_REQUIRED}")
set(_PGO_SAVED_CMAKE_CXX_EXTENSIONS "${CMAKE_CXX_EXTENSIONS}")

set(FLOAT_TETWILD_ENABLE_TBB ON CACHE BOOL "" FORCE)
set(FLOAT_TETWILD_USE_FLOAT OFF CACHE BOOL "" FORCE)
set(FLOAT_TETWILD_WITH_SANITIZERS OFF CACHE BOOL "" FORCE)
set(FLOAT_TETWILD_WITH_EXACT_ENVELOPE OFF CACHE BOOL "" FORCE)

function(_libpgo_patch_ftetwild_geogram target_file)
  file(READ "${target_file}" _libpgo_geogram_cmake)
  set(_libpgo_patched_uninstall [=[if(NOT TARGET uninstall)
add_custom_target(uninstall
COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()]=])
  if("${_libpgo_geogram_cmake}" MATCHES "if\\(NOT TARGET uninstall\\)[\r\n]+add_custom_target\\(uninstall")
    return()
  endif()

  set(_libpgo_unpatched_uninstall [=[add_custom_target(uninstall
COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)]=])
  string(FIND "${_libpgo_geogram_cmake}" "${_libpgo_unpatched_uninstall}" _libpgo_uninstall_index)
  if(_libpgo_uninstall_index EQUAL -1)
    message(FATAL_ERROR "Failed to patch ${target_file}: Geogram uninstall target snippet was not found.")
  endif()

  string(REPLACE
    "${_libpgo_unpatched_uninstall}"
    "${_libpgo_patched_uninstall}"
    _libpgo_geogram_cmake
    "${_libpgo_geogram_cmake}")
  file(WRITE "${target_file}" "${_libpgo_geogram_cmake}")
endfunction()

function(_libpgo_prepare_ftetwild_geogram)
  if(TARGET geogram)
    if(NOT TARGET geogram::geogram)
      add_library(geogram::geogram ALIAS geogram)
    endif()
    return()
  endif()

  if(MSVC)
    set(GEO_PLATFORM "Win-vs-generic")
  elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(GEO_PLATFORM "Darwin-clang")
  else()
    set(GEO_PLATFORM "Linux64-gcc")
  endif()

  set(VORPALINE_PLATFORM ${GEO_PLATFORM} CACHE STRING "" FORCE)
  set(GEOGRAM_BUILD_SHARED OFF CACHE BOOL "" FORCE)
  set(GEOGRAM_BUILD_STATIC ON CACHE BOOL "" FORCE)
  set(GEOGRAM_SUB_BUILD ON CACHE BOOL "Building as subproject" FORCE)
  set(GEOGRAM_LIB_ONLY ON CACHE BOOL "Build geogram lib only" FORCE)
  set(GEOGRAM_WITH_GRAPHICS OFF CACHE BOOL "Disable graphics" FORCE)
  set(GEOGRAM_WITH_LUA OFF CACHE BOOL "Disable LUA" FORCE)
  set(GEOGRAM_WITH_EXPLORAGRAM OFF CACHE BOOL "Disable exploragram" FORCE)
  set(GEOGRAM_WITH_LEGACY_NUMERICS OFF CACHE BOOL "Disable legacy numerics" FORCE)
  set(GEOGRAM_WITH_TRIANGLE OFF CACHE BOOL "Disable triangle" FORCE)

  FetchContent_Declare(
    geogram
    GIT_REPOSITORY https://github.com/BrunoLevy/geogram
    GIT_TAG v1.9.6
  )

  FetchContent_GetProperties(geogram)
  if(NOT geogram_POPULATED)
    FetchContent_Populate(geogram)
  endif()

  _libpgo_patch_ftetwild_geogram("${geogram_SOURCE_DIR}/CMakeLists.txt")
  add_subdirectory(${geogram_SOURCE_DIR} ${geogram_BINARY_DIR} EXCLUDE_FROM_ALL)

  if(TARGET geogram AND NOT TARGET geogram::geogram)
    add_library(geogram::geogram ALIAS geogram)
  endif()
endfunction()

# fTetWild fetches Geogram itself; preloading a patched target avoids duplicate
# global targets such as Eigen's and Geogram's "uninstall".
_libpgo_prepare_ftetwild_geogram()

FetchContent_Declare(
  ftetwild
  GIT_REPOSITORY https://github.com/wildmeshing/fTetWild.git
  GIT_TAG d7d99bb4387a07895b9adce058dc7305f6b6e5ab
)

FetchContent_MakeAvailable(ftetwild)

set(CMAKE_CXX_STANDARD "${_PGO_SAVED_CMAKE_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD_REQUIRED "${_PGO_SAVED_CMAKE_CXX_STANDARD_REQUIRED}")
set(CMAKE_CXX_EXTENSIONS "${_PGO_SAVED_CMAKE_CXX_EXTENSIONS}")

if(NOT TARGET FloatTetwild)
  message(FATAL_ERROR "PGO_TET_MESHER_USE_TET_WILD=ON requires fTetWild target FloatTetwild, but it was not created.")
endif()

set(_PGO_FTETWILD_COMPAT_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/pgo_ftetwild_compat/include")
file(MAKE_DIRECTORY "${_PGO_FTETWILD_COMPAT_INCLUDE_DIR}/igl/predicates")
file(WRITE "${_PGO_FTETWILD_COMPAT_INCLUDE_DIR}/igl/predicates/predicates.h"
  "#pragma once\n"
  "#include <igl/Orientation.h>\n"
  "namespace igl { namespace predicates { using Orientation = ::igl::Orientation; } }\n"
  "#include <igl/predicates/exactinit.h>\n"
  "#include <igl/predicates/orient2d.h>\n"
  "#include <igl/predicates/orient3d.h>\n")
target_include_directories(FloatTetwild PRIVATE "${_PGO_FTETWILD_COMPAT_INCLUDE_DIR}")

message(STATUS "Done.")
