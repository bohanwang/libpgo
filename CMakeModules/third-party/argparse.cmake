if(TARGET argparse::argparse)
  return()
endif()

message(STATUS "Loading argparse...")

set(ARGPARSE_INSTALL OFF CACHE BOOL "Include an install target" FORCE)
set(ARGPARSE_BUILD_TESTS OFF CACHE BOOL "Build tests" FORCE)
set(ARGPARSE_BUILD_SAMPLES OFF CACHE BOOL "Build samples" FORCE)

include(FetchContent)
FetchContent_Declare(
  argparse
  URL https://github.com/p-ranav/argparse/archive/refs/tags/v3.1.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES argparse
)

FetchContent_MakeAvailable(argparse)

message(STATUS "Done.")


