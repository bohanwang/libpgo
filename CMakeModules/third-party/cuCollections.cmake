if(TARGET cuco)
  return()
endif()

message(STATUS "Loading cuco...")

pgo_dep_option(BUILD_TESTS BOOL OFF "Configure CMake to build tests")
pgo_dep_option(BUILD_BENCHMARKS BOOL OFF "Configure CMake to build (google) benchmarks")
pgo_dep_option(BUILD_EXAMPLES BOOL OFF "Configure CMake to build examples")
pgo_dep_option(BUILD_CUCO_TESTS BOOL OFF "Configure CMake to build cuco tests")

include(FetchContent)
FetchContent_Declare(
  cuco
  GIT_REPOSITORY https://github.com/NVIDIA/cuCollections.git
  GIT_TAG dev
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(cuco)

message(STATUS "Done.")
