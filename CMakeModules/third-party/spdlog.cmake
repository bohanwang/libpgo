if(TARGET spdlog::spdlog_header_only)
  return()
endif()

message(STATUS "Loading spdlog...")

include(FetchContent)
FetchContent_Declare(
  spdlog
  URL https://github.com/gabime/spdlog/archive/refs/tags/v1.15.2.zip
  URL_HASH SHA256=d91ab0e16964cedb826e65ba1bed5ed4851d15c7b9453609a52056a94068c020
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(spdlog)

message(STATUS "Done.")
