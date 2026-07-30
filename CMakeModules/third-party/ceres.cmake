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
# libpgo links SuiteSparse directly. Letting Ceres discover the same in-tree
# targets makes its legacy find module try to mutate alias targets.
set(SUITESPARSE OFF CACHE BOOL "Use SuiteSparse in Ceres" FORCE)

include(FetchContent)
FetchContent_Declare(
  ceres
  URL https://github.com/ceres-solver/ceres-solver/archive/refs/tags/2.2.0.tar.gz
  URL_HASH SHA256=12efacfadbfdc1bbfa203c236e96f4d3c210bed96994288b3ff0c8e7c6f350d4
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(ceres)

message(STATUS "Done.")
