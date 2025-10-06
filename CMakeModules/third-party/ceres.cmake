if(TARGET Ceres::Ceres)
  return()
endif()

message(STATUS "Loading Ceres...")

set(USE_CUDA OFF CACHE BOOL "use cuda" FORCE)
set(MINIGLOG ON CACHE BOOL "use mini glog" FORCE)
set(GFLAGS OFF CACHE BOOL "use gflags" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Enable tests" FORCE)
set(BUILD_DOCUMENTATION OFF CACHE BOOL "Build User's Guide (html)" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "Build examples" FORCE)
set(BUILD_BENCHMARKS OFF CACHE BOOL "Build Ceres benchmarking suite" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Ceres as a shared library." FORCE)
set(PROVIDE_UNINSTALL_TARGET OFF CACHE BOOL "Add a custom target to ease removal of installed targets" FORCE)
set(LAPACK OFF CACHE BOOL "Use LAPACK" FORCE)

include(FetchContent)
FetchContent_Declare(
  ceres
  URL http://ceres-solver.org/ceres-solver-2.2.0.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES Ceres
)

FetchContent_MakeAvailable(ceres)

message(STATUS "Done.")
