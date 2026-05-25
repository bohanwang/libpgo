if(TARGET spdlog::spdlog_header_only)
  return()
endif()

message(STATUS "Loading spdlog...")

include(FetchContent)
pgo_dep_option(SPDLOG_FMT_EXTERNAL_HO BOOL ON "Use the existing header-only fmt target in spdlog")
FetchContent_Declare(
  spdlog
  URL https://github.com/gabime/spdlog/archive/refs/tags/v1.15.2.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(spdlog)

message(STATUS "Done.")
