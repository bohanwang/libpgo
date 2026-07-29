if(TARGET cuco)
  return()
endif()

message(STATUS "Loading cuco...")

set(BUILD_TESTS OFF CACHE BOOL "Configure CMake to build tests" FORCE)
set(BUILD_BENCHMARKS OFF CACHE BOOL "Configure CMake to build (google) benchmarks" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "Configure CMake to build examples" FORCE)
set(BUILD_CUCO_TESTS OFF CACHE BOOL "Configure CMake to build cuco tests" FORCE)

include(FetchContent)
FetchContent_Declare(
  cuco
  GIT_REPOSITORY https://github.com/NVIDIA/cuCollections.git
  GIT_TAG 1a1e640179e85139765a27d2d376e02628b2ccbc
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES cuco
)

FetchContent_MakeAvailable(cuco)

message(STATUS "Done.")
