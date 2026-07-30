option(
  PGO_STAGE_WINDOWS_RUNTIME
  "Stage Windows runtime DLLs beside source-build executables and pypgo"
  ON
)

function(pgo_stage_windows_runtime_for_module target)
  if(NOT TARGET pgo_windows_runtime)
    return()
  endif()

  add_dependencies(${target} pgo_windows_runtime)
  add_custom_command(
    TARGET ${target}
    POST_BUILD
    COMMAND
      ${CMAKE_COMMAND} -E copy_directory_if_different
      "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}"
      "$<TARGET_FILE_DIR:${target}>"
    COMMENT "Staging Windows runtime DLLs beside $<TARGET_FILE_NAME:${target}>"
    VERBATIM
  )
endfunction()

function(_pgo_collect_project_runtime_targets directory)
  get_property(_pgo_targets DIRECTORY "${directory}" PROPERTY BUILDSYSTEM_TARGETS)
  foreach(_pgo_target IN LISTS _pgo_targets)
    get_target_property(_pgo_target_imported "${_pgo_target}" IMPORTED)
    if(_pgo_target_imported)
      continue()
    endif()

    get_target_property(_pgo_target_source_dir "${_pgo_target}" SOURCE_DIR)
    set(_pgo_source_root "${CMAKE_SOURCE_DIR}")
    cmake_path(
      IS_PREFIX
      _pgo_source_root
      "${_pgo_target_source_dir}"
      NORMALIZE
      _pgo_is_project_target
    )
    if(NOT _pgo_is_project_target)
      continue()
    endif()

    get_target_property(_pgo_target_type "${_pgo_target}" TYPE)
    if(
      NOT _pgo_target_type STREQUAL "EXECUTABLE"
      AND NOT _pgo_target_type STREQUAL "SHARED_LIBRARY"
    )
      continue()
    endif()

    set_target_properties(
      ${_pgo_target}
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
          "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG
          "${CMAKE_BINARY_DIR}/bin/Debug"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE
          "${CMAKE_BINARY_DIR}/bin/Release"
        RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO
          "${CMAKE_BINARY_DIR}/bin/RelWithDebInfo"
        RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL
          "${CMAKE_BINARY_DIR}/bin/MinSizeRel"
    )
    add_dependencies(${_pgo_target} pgo_windows_runtime)
  endforeach()

  get_property(
    _pgo_subdirectories
    DIRECTORY "${directory}"
    PROPERTY SUBDIRECTORIES
  )
  foreach(_pgo_subdirectory IN LISTS _pgo_subdirectories)
    _pgo_collect_project_runtime_targets("${_pgo_subdirectory}")
  endforeach()
endfunction()

function(pgo_finalize_windows_runtime_targets)
  if(NOT TARGET pgo_windows_runtime)
    return()
  endif()

  _pgo_collect_project_runtime_targets("${CMAKE_SOURCE_DIR}")
endfunction()

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

set(_pgo_windows_runtime_dlls
  "${_pgo_windows_python_bin}/tbb12.dll"
  "${CMAKE_SOURCE_DIR}/third-party/gmp-msvc/release/gmp-10.dll"
  "${CMAKE_SOURCE_DIR}/third-party/gmp-msvc/release/gmpxx-4.dll"
  "${CMAKE_SOURCE_DIR}/third-party/mpfr-msvc/release/mpfr-6.dll"
)

if(PGO_USE_MKL)
  foreach(_pgo_mkl_runtime_name
      core
      tbb_thread
      def
      mc3
      avx2
      avx512
      vml_def
      vml_cmpt
      vml_mc3
      vml_avx2
      vml_avx512)
    list(APPEND
      _pgo_windows_runtime_dlls
      "${_pgo_windows_python_bin}/mkl_${_pgo_mkl_runtime_name}.2.dll"
    )
  endforeach()
endif()

foreach(_pgo_windows_required_dll IN LISTS _pgo_windows_runtime_dlls)
  if(NOT EXISTS "${_pgo_windows_required_dll}")
    message(FATAL_ERROR
      "Required Windows runtime DLL is missing: "
      "${_pgo_windows_required_dll}")
  endif()
endforeach()

list(REMOVE_DUPLICATES _pgo_windows_runtime_dlls)
list(SORT _pgo_windows_runtime_dlls)

set(
  PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY
  "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
  CACHE INTERNAL
  "Shared Windows source-build runtime directory"
  FORCE
)

set(_pgo_windows_runtime_manifest "")
foreach(_pgo_windows_runtime_dll IN LISTS _pgo_windows_runtime_dlls)
  string(APPEND
    _pgo_windows_runtime_manifest
    "${_pgo_windows_runtime_dll}\n")
endforeach()
file(
  WRITE
  "${CMAKE_BINARY_DIR}/windows-runtime-manifest.txt"
  "${_pgo_windows_runtime_manifest}"
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

message(STATUS
  "Windows runtime staging: ${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}")
