if(TARGET libcmaes::cmaes)
  return()
endif()

message(STATUS "Loading libcmaes...")

set(LIBCMAES_BUILD_SHARED_LIBS OFF CACHE BOOL "Build libcmaes as a shared library" FORCE)
set(LIBCMAES_BUILD_PYTHON OFF CACHE BOOL "build python bindings" FORCE)
set(LIBCMAES_BUILD_TESTS OFF CACHE BOOL "Build tests" FORCE)
set(LIBCMAES_BUILD_EXAMPLES OFF CACHE BOOL "Build samples" FORCE)
set(LIBCMAES_USE_OPENMP OFF CACHE BOOL "Use OpenMP for multithreading" FORCE)
set(LIBCMAES_ENABLE_SURROG ON CACHE BOOL "support for surrogates"  FORCE)

include(FetchContent)
FetchContent_Declare(
  libcmaes
  GIT_REPOSITORY https://github.com/bohanwang/libcmaes.git
  GIT_TAG 4598001f8174a0f2b90fc8983f819c8d047af26f
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES libcmaes
)

FetchContent_MakeAvailable(libcmaes)

message(STATUS "Done.")


