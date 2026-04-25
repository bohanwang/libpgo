message(STATUS "Loading fTetWild...")

include(FetchContent)

set(_PGO_SAVED_CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD}")
set(_PGO_SAVED_CMAKE_CXX_STANDARD_REQUIRED "${CMAKE_CXX_STANDARD_REQUIRED}")
set(_PGO_SAVED_CMAKE_CXX_EXTENSIONS "${CMAKE_CXX_EXTENSIONS}")

set(FLOAT_TETWILD_ENABLE_TBB ON CACHE BOOL "" FORCE)
set(FLOAT_TETWILD_USE_FLOAT OFF CACHE BOOL "" FORCE)
set(FLOAT_TETWILD_WITH_SANITIZERS OFF CACHE BOOL "" FORCE)
set(FLOAT_TETWILD_WITH_EXACT_ENVELOPE OFF CACHE BOOL "" FORCE)

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
