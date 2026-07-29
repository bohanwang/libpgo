if(TARGET SuiteSparse::SuiteSparse_config)
  return()
endif()

message(STATUS "Loading SuiteSparse...")

set(SUITESPARSE_USE_CUDA OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_DEMOS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_ENABLE_PROJECTS "suitesparse_config;amd;camd;ccolamd;colamd;cholmod;cxsparse;klu;umfpack;spqr;" CACHE STRING "" FORCE)
set(SUITESPARSE_USE_FORTRAN OFF CACHE BOOL "" FORCE)
set(SUITESPARSE_USE_OPENMP OFF CACHE BOOL "" FORCE)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(APPLE AND NOT PGO_USE_MKL)
  get_target_property(PGO_BLAS_LINK_LIBRARIES BLAS::BLAS INTERFACE_LINK_LIBRARIES)
  get_target_property(PGO_LAPACK_LINK_LIBRARIES LAPACK::LAPACK INTERFACE_LINK_LIBRARIES)
  string(JOIN ";" PGO_BLAS_LINK_INTERFACE ${PGO_BLAS_LINK_LIBRARIES})
  string(JOIN ";" PGO_LAPACK_LINK_INTERFACE ${PGO_LAPACK_LINK_LIBRARIES})
  string(JOIN ";" PGO_LAPACK_PROVIDER
    ${PGO_LAPACK_LINK_LIBRARIES}
    ${LAPACK_LIBRARIES})
  if(NOT PGO_BLAS_LINK_INTERFACE MATCHES "Accelerate"
      OR NOT PGO_LAPACK_PROVIDER MATCHES "Accelerate")
    message(FATAL_ERROR
      "macOS release builds require system Accelerate for BLAS and LAPACK. "
      "BLAS=${PGO_BLAS_LINK_INTERFACE}; LAPACK=${PGO_LAPACK_PROVIDER}")
  endif()
  message(STATUS "macOS BLAS provider: ${PGO_BLAS_LINK_INTERFACE}")
  message(STATUS "macOS LAPACK provider: ${PGO_LAPACK_PROVIDER}")
endif()

include(FetchContent)
FetchContent_Declare(
  suitesparse
  URL https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v7.10.3.zip
  URL_HASH SHA256=d4600765554133fb3c0a830ace87ff89225a01ffe02f4012ecfdcf7c778dbcc0
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES SuiteSparse COMPONENTS SuiteSparse_config cholmod spqr umfpack
)

FetchContent_MakeAvailable(suitesparse)

message(STATUS "Done.")
