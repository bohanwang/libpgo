if(TARGET fmt::fmt-header-only)
  return()
endif()

message(STATUS "Loading fmtlib...")

include(FetchContent)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/11.1.4.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES fmt
)

FetchContent_MakeAvailable(fmt)

message(STATUS "Done.")