if(TARGET libcmaes::cmaes)
  return()
endif()

message(STATUS "Loading libcmaes...")

pgo_dep_option(LIBCMAES_BUILD_SHARED_LIBS BOOL OFF "Build libcmaes as a shared library")
pgo_dep_option(LIBCMAES_BUILD_PYTHON BOOL OFF "build python bindings")
pgo_dep_option(LIBCMAES_BUILD_TESTS BOOL OFF "Build tests")
pgo_dep_option(LIBCMAES_BUILD_EXAMPLES BOOL OFF "Build samples")
pgo_dep_option(LIBCMAES_USE_OPENMP BOOL OFF "Use OpenMP for multithreading")
pgo_dep_option(LIBCMAES_ENABLE_SURROG BOOL ON "support for surrogates")

include(FetchContent)
FetchContent_Declare(
  libcmaes
  GIT_REPOSITORY https://github.com/bohanwang/libcmaes.git
  GIT_TAG 4598001f8174a0f2b90fc8983f819c8d047af26f
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(libcmaes)

message(STATUS "Done.")

