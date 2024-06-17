if(TARGET pybind11::pybind11)
  return()
endif()

message(STATUS "Loading pybind11...")

set(PYBIND11_TEST OFF)

include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.12.0.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES pybind11
)

FetchContent_MakeAvailable(pybind11)

message(STATUS "Done.")