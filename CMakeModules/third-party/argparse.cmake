if(TARGET argparse::argparse)
  return()
endif()

message(STATUS "Loading argparse...")

pgo_dep_option(ARGPARSE_INSTALL BOOL OFF "Include an install target")
pgo_dep_option(ARGPARSE_BUILD_TESTS BOOL OFF "Build tests")
pgo_dep_option(ARGPARSE_BUILD_SAMPLES BOOL OFF "Build samples")

include(FetchContent)
FetchContent_Declare(
  argparse
  URL https://github.com/p-ranav/argparse/archive/refs/tags/v3.1.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(argparse)

message(STATUS "Done.")

