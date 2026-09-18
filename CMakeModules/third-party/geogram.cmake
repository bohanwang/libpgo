if(TARGET geogram::geogram)
  return()
endif()

message(STATUS "Loading geogram...")

set(GEOGRAM_SUB_BUILD ON CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_HLBFGS ON CACHE BOOL "Non-linear solver (Yang Liu's HLBFGS)" FORCE)

include(FetchContent)
FetchContent_Declare(
  geogram
  URL https://github.com/BrunoLevy/geogram/releases/download/v1.9.0/geogram_1.9.0.zip
  URL_HASH SHA256=f2b51adf05fc8599893032c79866b4f2ff29326f810dcde649adf205896b76ab
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES geogram
)

FetchContent_GetProperties(geogram)
if(NOT geogram_POPULATED)
  FetchContent_Populate(geogram)
endif()

function(_libpgo_replace_in_file target_file old_text new_text)
  file(READ "${target_file}" _libpgo_file_contents)
  string(FIND "${_libpgo_file_contents}" "${new_text}" _libpgo_already_patched_index)
  if(NOT _libpgo_already_patched_index EQUAL -1)
    return()
  endif()
  string(FIND "${_libpgo_file_contents}" "${old_text}" _libpgo_match_index)
  if(_libpgo_match_index EQUAL -1)
    message(FATAL_ERROR "Failed to patch ${target_file}: expected snippet was not found.")
  endif()
  string(REPLACE "${old_text}" "${new_text}" _libpgo_file_contents "${_libpgo_file_contents}")
  file(WRITE "${target_file}" "${_libpgo_file_contents}")
endfunction()

set(MODIFIED_FILE "${CMAKE_SOURCE_DIR}/CMakeModules/patches/geogram.cmake")
set(TARGET_FILE "${geogram_SOURCE_DIR}/CMakeLists.txt")

# Read in the content
file(READ "${MODIFIED_FILE}" content)

# Write the modified content back to the file
file(WRITE "${TARGET_FILE}" "${content}")

set(POISSON_RECON_DIR "${geogram_SOURCE_DIR}/src/lib/geogram/third_party/PoissonRecon")

_libpgo_replace_in_file(
  "${POISSON_RECON_DIR}/SparseMatrix.inl"
  [=[void SparseMatrix<T>::SetZero()
{
        Resize(this->m_N, this->m_M);
}]=]
  [=[void SparseMatrix<T>::SetZero()
{
        for( int i=0 ; i<rows ; i++ ) for( int ii=0 ; ii<rowSizes[i] ; ii++ ) m_ppElements[i][ii].Value = T(0);
}]=]
)

_libpgo_replace_in_file(
  "${POISSON_RECON_DIR}/SparseMatrix.inl"
  [=[for( int i=0 ; i<rows ; i++ ) for( int ii=0 ; ii<rowSizes[i] ; i++ ) m_ppElements[i][ii].Value *= V;]=]
  [=[for( int i=0 ; i<rows ; i++ ) for( int ii=0 ; ii<rowSizes[i] ; ii++ ) m_ppElements[i][ii].Value *= V;]=]
)

_libpgo_replace_in_file(
  "${POISSON_RECON_DIR}/PlyVertexMini.h"
  [=[        PlyValueVertex operator - ( PlyValueVertex p ) const { return PlyValueVertex( point-p.value , value-p.value ); }]=]
  [=[        PlyValueVertex operator - ( PlyValueVertex p ) const { return PlyValueVertex( point-p.point , value-p.value ); }]=]
)

_libpgo_replace_in_file(
  "${POISSON_RECON_DIR}/PlyVertexMini.h"
  [=[        PlyOrientedVertex operator - ( PlyOrientedVertex p ) const { return PlyOrientedVertex( point-p.value , normal-p.normal ); }]=]
  [=[        PlyOrientedVertex operator - ( PlyOrientedVertex p ) const { return PlyOrientedVertex( point-p.point , normal-p.normal ); }]=]
)

_libpgo_replace_in_file(
  "${POISSON_RECON_DIR}/PlyVertexMini.h"
  [=[                _PlyColorVertex operator - ( _PlyColorVertex p ) const { return _PlyColorVertex( point-p.value , color-p.color ); }]=]
  [=[                _PlyColorVertex operator - ( _PlyColorVertex p ) const { return _PlyColorVertex( point-p.point , color-p.color ); }]=]
)

_libpgo_replace_in_file(
  "${POISSON_RECON_DIR}/PlyVertexMini.h"
  [=[                _PlyColorAndValueVertex operator - ( _PlyColorAndValueVertex p ) const { return _PlyColorAndValueVertex( point-p.value , color-p.color , value+p.value ); }]=]
  [=[                _PlyColorAndValueVertex operator - ( _PlyColorAndValueVertex p ) const { return _PlyColorAndValueVertex( point-p.point , color-p.color , value-p.value ); }]=]
)

add_subdirectory(${geogram_SOURCE_DIR} ${geogram_BINARY_DIR} EXCLUDE_FROM_ALL)

if(TARGET OpenMP::OpenMP_CXX)
  # HLBFGS is compiled as an OBJECT library and directly includes omp.h.
  target_link_libraries(geogram_third_party PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(geogram_third_party PRIVATE USE_OPENMP)

  # Geogram itself uses the legacy plain target_link_libraries signature.
  # Keep the OpenMP runtime on the final Geogram link interface so consumers
  # also receive it when Geogram is built from object/static sources.
  target_link_libraries(geogram OpenMP::OpenMP_CXX)
endif()

set_target_properties(geogram PROPERTIES CXX_STANDARD 14)

message(STATUS "Done.")
