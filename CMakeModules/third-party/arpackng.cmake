if(TARGET ARPACK::ARPACK)
  return()
endif()

message(STATUS "Loading arpack...")

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries instead of static libraries" FORCE)
set(MPI OFF CACHE BOOL "Enable parallel support" FORCE)
set(ICB OFF CACHE BOOL "Enable support for *[ae]upd_c with ISO_C_BINDING" FORCE)
set(EIGEN OFF CACHE BOOL "Enable support for eigenvalue-problems solver based on ICB and eigen" FORCE)
set(PYTHON3 OFF CACHE BOOL "Enable python3 support" FORCE)
set(EXAMPLES OFF CACHE BOOL "Compile ARPACK examples" FORCE)
set(TESTS OFF CACHE BOOL "Compile ARPACK tests" FORCE)


include(FetchContent)
FetchContent_Declare(
  arpackng
  GIT_REPOSITORY https://github.com/opencollab/arpack-ng.git
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES arpackng
)

FetchContent_MakeAvailable(arpackng)

FetchContent_GetProperties(arpackng)
if (NOT arpackng_POPULATED)
  message(FATAL_ERROR "Failed to download arpack")
else()
  set(PGO_HAS_ARPACK ON CACHE BOOL "ARPACK is available" FORCE)
endif()

message(STATUS "Done.")



