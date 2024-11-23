if(TARGET spdlog::spdlog_header_only)
  return()
endif()

message(STATUS "Loading spdlog...")

include(FetchContent)
FetchContent_Declare(
  spdlog
  URL https://github.com/gabime/spdlog/archive/refs/tags/v1.13.0.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES spdlog
)

FetchContent_MakeAvailable(spdlog)

message(STATUS "Done.")