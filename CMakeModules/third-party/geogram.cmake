if(TARGET geogram::geogram)
  return()
endif()

message(STATUS "Loading geogram...")

pgo_dep_option(GEOGRAM_SUB_BUILD BOOL ON "Building as subproject")
pgo_dep_option(GEOGRAM_LIB_ONLY BOOL ON "Build geogram lib only")
pgo_dep_option(GEOGRAM_WITH_GRAPHICS BOOL OFF "Disable graphics")
pgo_dep_option(GEOGRAM_WITH_HLBFGS BOOL ON "Non-linear solver (Yang Liu's HLBFGS)")
pgo_dep_option(GEOGRAM_WITH_LUA BOOL OFF "Disable LUA")
pgo_dep_option(GEOGRAM_WITH_EXPLORAGRAM BOOL OFF "Disable exploragram")
pgo_dep_option(GEOGRAM_WITH_LEGACY_NUMERICS BOOL OFF "Disable legacy numerics")
pgo_dep_option(GEOGRAM_WITH_TRIANGLE BOOL OFF "Disable triangle")

include(FetchContent)
FetchContent_Declare(
  geogram
  URL https://github.com/BrunoLevy/geogram/releases/download/v1.9.0/geogram_1.9.0.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_populate_compat(geogram "geogram source tree is patched before add_subdirectory")

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

file(READ "${MODIFIED_FILE}" content)
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

set_target_properties(geogram PROPERTIES CXX_STANDARD 14)

message(STATUS "Done.")
