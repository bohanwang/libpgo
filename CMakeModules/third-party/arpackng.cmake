if(TARGET ARPACK::ARPACK)
  return()
endif()

message(STATUS "Loading arpack...")

pgo_dep_option(BUILD_SHARED_LIBS BOOL OFF "Build shared libraries instead of static libraries")
pgo_dep_option(MPI BOOL OFF "Enable parallel support")
pgo_dep_option(ICB BOOL OFF "Enable support for *[ae]upd_c with ISO_C_BINDING")
pgo_dep_option(EIGEN BOOL OFF "Enable support for eigenvalue-problems solver based on ICB and eigen")
pgo_dep_option(PYTHON3 BOOL OFF "Enable python3 support")
pgo_dep_option(EXAMPLES BOOL OFF "Compile ARPACK examples")
pgo_dep_option(TESTS BOOL OFF "Compile ARPACK tests")


include(FetchContent)
FetchContent_Declare(
  arpackng
  GIT_REPOSITORY https://github.com/opencollab/arpack-ng.git
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(arpackng)

FetchContent_GetProperties(arpackng)
if (NOT arpackng_POPULATED)
  message(FATAL_ERROR "Failed to download arpack")
else()
  pgo_dep_option(PGO_HAS_ARPACK BOOL ON "ARPACK is available")
endif()

message(STATUS "Done.")


