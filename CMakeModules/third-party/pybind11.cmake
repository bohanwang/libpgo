if(TARGET pybind11::pybind11)
  return()
endif()

message(STATUS "Loading pybind11...")

set(PYBIND11_TEST OFF)

include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.12.0.tar.gz
  URL_HASH SHA256=bf8f242abd1abcd375d516a7067490fb71abd79519a282d22b6e4d19282185a7
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(pybind11)

message(STATUS "Done.")
