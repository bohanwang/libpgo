if(TARGET Ceres::Ceres)
  return()
endif()

message(STATUS "Loading Ceres...")

pgo_dep_option(USE_CUDA BOOL OFF "use cuda")
pgo_dep_option(MINIGLOG BOOL ON "use mini glog")
pgo_dep_option(GFLAGS BOOL OFF "use gflags")
pgo_dep_option(BUILD_TESTING BOOL OFF "Enable tests")
pgo_dep_option(BUILD_DOCUMENTATION BOOL OFF "Build User's Guide (html)")
pgo_dep_option(BUILD_EXAMPLES BOOL OFF "Build examples")
pgo_dep_option(BUILD_BENCHMARKS BOOL OFF "Build Ceres benchmarking suite")
pgo_dep_option(BUILD_SHARED_LIBS BOOL OFF "Build Ceres as a shared library.")
pgo_dep_option(PROVIDE_UNINSTALL_TARGET BOOL OFF "Add a custom target to ease removal of installed targets")
pgo_dep_option(LAPACK BOOL OFF "Use LAPACK")

include(FetchContent)
FetchContent_Declare(
  ceres
  URL http://ceres-solver.org/ceres-solver-2.2.0.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(ceres)

message(STATUS "Done.")
