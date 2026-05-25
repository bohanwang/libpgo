if(TARGET SuiteSparse::SuiteSparse_config)
  return()
endif()

message(STATUS "Loading SuiteSparse...")

pgo_dep_option(SUITESPARSE_USE_CUDA BOOL OFF "SuiteSparse CUDA support")
pgo_dep_option(SUITESPARSE_DEMOS BOOL OFF "SuiteSparse demos")
pgo_dep_option(BUILD_SHARED_LIBS BOOL OFF "Build shared libraries")
pgo_dep_option(SUITESPARSE_ENABLE_PROJECTS STRING "suitesparse_config;amd;camd;ccolamd;colamd;cholmod;cxsparse;klu;umfpack;spqr;" "SuiteSparse projects to build")
pgo_dep_option(SUITESPARSE_USE_FORTRAN BOOL OFF "SuiteSparse Fortran support")
pgo_dep_option(SUITESPARSE_USE_OPENMP BOOL OFF "SuiteSparse OpenMP support")
pgo_dep_option(BLA_VENDOR STRING "" "BLAS vendor")

include(FetchContent)
FetchContent_Declare(
  suitesparse
  URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.10.3.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(suitesparse)

foreach(_pgo_suitesparse_component IN ITEMS AMD CAMD CCOLAMD CHOLMOD COLAMD SPQR)
  if(TARGET ${_pgo_suitesparse_component}_static AND NOT TARGET SuiteSparse::${_pgo_suitesparse_component})
    add_library(SuiteSparse::${_pgo_suitesparse_component} ALIAS ${_pgo_suitesparse_component}_static)
  endif()
endforeach()

if(TARGET SuiteSparseConfig_static)
  if(NOT TARGET SuiteSparse::Config)
    add_library(SuiteSparse::Config ALIAS SuiteSparseConfig_static)
  endif()
  if(NOT TARGET SuiteSparse::SuiteSparseConfig)
    add_library(SuiteSparse::SuiteSparseConfig ALIAS SuiteSparseConfig_static)
  endif()
endif()

# Ceres calls find_package(SuiteSparse) from inside its own FetchContent build.
# Point that lookup back at the fetched SuiteSparse targets instead of letting it
# scan system locations.
file(WRITE "${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/SuiteSparseConfig.cmake"
  [=[
set(SuiteSparse_FOUND TRUE)
set(SuiteSparse_VERSION "7.10.3")
set(SuiteSparse_VERSION_MAJOR 7)
set(SuiteSparse_VERSION_MINOR 10)
set(SuiteSparse_VERSION_PATCH 3)
foreach(_pgo_component IN ITEMS AMD CAMD CCOLAMD CHOLMOD COLAMD SPQR Config)
  if(TARGET SuiteSparse::${_pgo_component})
    set(SuiteSparse_${_pgo_component}_FOUND TRUE)
  else()
    set(SuiteSparse_${_pgo_component}_FOUND FALSE)
    set(SuiteSparse_FOUND FALSE)
  endif()
endforeach()
set(SuiteSparse_Partition_FOUND FALSE)
]=])
file(WRITE "${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/SuiteSparseConfigVersion.cmake"
  [=[
set(PACKAGE_VERSION "7.10.3")
if(PACKAGE_FIND_VERSION VERSION_LESS_EQUAL PACKAGE_VERSION)
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
  if(PACKAGE_FIND_VERSION VERSION_EQUAL PACKAGE_VERSION)
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
]=])

message(STATUS "Done.")
