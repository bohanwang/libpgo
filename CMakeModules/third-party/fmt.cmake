if(TARGET fmt::fmt-header-only)
  return()
endif()

message(STATUS "Loading fmtlib...")

include(FetchContent)
FetchContent_Declare(
  fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/11.1.4.zip
  URL_HASH SHA256=7e85cbf6125a76daa0f83cd9240eff863d988aca68cd5f66c01ff7b59fa886b6
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(fmt)

message(STATUS "Done.")
