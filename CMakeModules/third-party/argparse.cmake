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
  URL_HASH SHA256=3e5a59ab7688dcd1f918bc92051a10564113d4f36c3bbed3ef596c25e519a062
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(argparse)

message(STATUS "Done.")
