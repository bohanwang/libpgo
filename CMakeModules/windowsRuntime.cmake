option(
  PGO_STAGE_WINDOWS_RUNTIME
  "Stage Windows runtime DLLs beside source-build executables and pypgo"
  ON
)

if(NOT WIN32 OR NOT PGO_STAGE_WINDOWS_RUNTIME)
  return()
endif()

if(NOT PGO_PYTHON_DEPENDENCY_PREFIX)
  message(FATAL_ERROR
    "Windows runtime staging requires an active Python environment. "
    "Run `uv sync --locked` and configure through `uv run cmake`.")
endif()

set(
  _pgo_windows_python_bin
  "${PGO_PYTHON_DEPENDENCY_PREFIX}/Library/bin"
)
if(NOT IS_DIRECTORY "${_pgo_windows_python_bin}")
  message(FATAL_ERROR
    "Windows native dependency directory does not exist: "
    "${_pgo_windows_python_bin}")
endif()

file(GLOB _pgo_tbb_runtime_dlls
  LIST_DIRECTORIES FALSE
  "${_pgo_windows_python_bin}/tbb*.dll")
if(NOT _pgo_tbb_runtime_dlls)
  message(FATAL_ERROR
    "No TBB runtime DLLs found in ${_pgo_windows_python_bin}")
endif()

set(_pgo_windows_runtime_dlls ${_pgo_tbb_runtime_dlls})

if(PGO_USE_MKL)
  file(GLOB _pgo_mkl_runtime_dlls
    LIST_DIRECTORIES FALSE
    "${_pgo_windows_python_bin}/mkl_*.dll")
  if(NOT _pgo_mkl_runtime_dlls)
    message(FATAL_ERROR
      "No MKL runtime DLLs found in ${_pgo_windows_python_bin}")
  endif()
  list(APPEND _pgo_windows_runtime_dlls ${_pgo_mkl_runtime_dlls})
endif()

set(_pgo_project_runtime_dlls
  "${CMAKE_SOURCE_DIR}/third-party/gmp-msvc/release/gmp-10.dll"
  "${CMAKE_SOURCE_DIR}/third-party/gmp-msvc/release/gmpxx-4.dll"
  "${CMAKE_SOURCE_DIR}/third-party/mpfr-msvc/release/mpfr-6.dll"
)
foreach(_pgo_windows_required_dll IN LISTS _pgo_project_runtime_dlls)
  if(NOT EXISTS "${_pgo_windows_required_dll}")
    message(FATAL_ERROR
      "Project runtime DLL is missing: ${_pgo_windows_required_dll}")
  endif()
endforeach()
list(APPEND _pgo_windows_runtime_dlls ${_pgo_project_runtime_dlls})

list(REMOVE_DUPLICATES _pgo_windows_runtime_dlls)
list(SORT _pgo_windows_runtime_dlls)

set(
  PGO_WINDOWS_RUNTIME_DLLS
  "${_pgo_windows_runtime_dlls}"
  CACHE INTERNAL
  "Windows source-build runtime DLLs"
  FORCE
)

set(
  PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
  CACHE INTERNAL
  "Shared Windows source-build runtime directory"
  FORCE
)

add_custom_target(
  pgo_windows_runtime
  ALL
  COMMAND
    ${CMAKE_COMMAND} -E make_directory
    "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}"
  COMMAND
    ${CMAKE_COMMAND} -E copy_if_different
    ${_pgo_windows_runtime_dlls}
    "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}"
  COMMENT "Staging shared Windows runtime DLLs"
  COMMAND_EXPAND_LISTS
  VERBATIM
)

list(LENGTH _pgo_windows_runtime_dlls _pgo_windows_runtime_dll_count)
message(STATUS "Windows runtime staging: "
  "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY} "
  "(${_pgo_windows_runtime_dll_count} DLLs)")
