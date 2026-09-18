include_guard(GLOBAL)

function(pgo_configure_runtime_dependencies)
  if(NOT PGO_RUNTIME_LAYOUT MATCHES "^(SOURCE|WHEEL)$")
    message(FATAL_ERROR
      "PGO_RUNTIME_LAYOUT must be SOURCE or WHEEL, got: "
      "${PGO_RUNTIME_LAYOUT}")
  endif()

  if(PGO_RUNTIME_LAYOUT STREQUAL "WHEEL")
    if(NOT PGO_ENABLE_PYTHON)
      message(FATAL_ERROR
        "PGO_RUNTIME_LAYOUT=WHEEL requires PGO_ENABLE_PYTHON=ON")
    endif()
    message(STATUS
      "Runtime dependency layout: WHEEL (native dependencies are repaired "
      "by the platform wheel packager)")
    return()
  endif()

  message(STATUS "Runtime dependency layout: SOURCE")

  if(UNIX AND NOT APPLE AND DEFINED PGO_PYTHON_DEPENDENCY_PREFIX)
    set(_pgo_python_runtime_dir "${PGO_PYTHON_DEPENDENCY_PREFIX}/lib")
    if(IS_DIRECTORY "${_pgo_python_runtime_dir}")
      set(_pgo_build_rpath ${CMAKE_BUILD_RPATH})
      list(APPEND _pgo_build_rpath "${_pgo_python_runtime_dir}")
      list(REMOVE_DUPLICATES _pgo_build_rpath)
      set(CMAKE_BUILD_RPATH "${_pgo_build_rpath}" PARENT_SCOPE)
      message(STATUS
        "Using Python dependency runtime path: ${_pgo_python_runtime_dir}")
      unset(_pgo_build_rpath)
    endif()
    unset(_pgo_python_runtime_dir)
  endif()

  if(NOT WIN32)
    return()
  endif()

  if(NOT PGO_PYTHON_DEPENDENCY_PREFIX)
    message(FATAL_ERROR
      "Windows source runtime staging requires an active Python environment. "
      "Run `uv sync --locked` and configure through `uv run cmake`.")
  endif()

  set(_pgo_windows_python_bin
    "${PGO_PYTHON_DEPENDENCY_PREFIX}/Library/bin")
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
    "${CMAKE_SOURCE_DIR}/third-party/mpfr-msvc/release/mpfr-6.dll")
  foreach(_pgo_windows_required_dll IN LISTS _pgo_project_runtime_dlls)
    if(NOT EXISTS "${_pgo_windows_required_dll}")
      message(FATAL_ERROR
        "Project runtime DLL is missing: ${_pgo_windows_required_dll}")
    endif()
  endforeach()
  list(APPEND _pgo_windows_runtime_dlls ${_pgo_project_runtime_dlls})

  list(REMOVE_DUPLICATES _pgo_windows_runtime_dlls)
  list(SORT _pgo_windows_runtime_dlls)

  set(PGO_WINDOWS_RUNTIME_DLLS
    "${_pgo_windows_runtime_dlls}"
    CACHE INTERNAL "Windows source-build runtime DLLs" FORCE)
  set(PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/bin/$<CONFIG>"
    CACHE INTERNAL "Shared Windows source-build runtime directory" FORCE)

  add_custom_target(pgo_windows_runtime ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory
      "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      ${_pgo_windows_runtime_dlls}
      "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}"
    COMMENT "Staging shared Windows source-build runtime DLLs"
    COMMAND_EXPAND_LISTS
    VERBATIM)

  list(LENGTH _pgo_windows_runtime_dlls _pgo_windows_runtime_dll_count)
  message(STATUS "Windows source runtime staging: "
    "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY} "
    "(${_pgo_windows_runtime_dll_count} DLLs)")
endfunction()

function(pgo_attach_runtime_dependencies target)
  set(_pgo_options COPY_RUNTIME WHEEL_EXTENSION)
  cmake_parse_arguments(PGO_ATTACH "${_pgo_options}" "" "" ${ARGN})

  if(NOT TARGET ${target})
    message(FATAL_ERROR
      "pgo_attach_runtime_dependencies called for unknown target: ${target}")
  endif()

  if(PGO_RUNTIME_LAYOUT STREQUAL "WHEEL")
    if(NOT PGO_ATTACH_WHEEL_EXTENSION)
      return()
    endif()

    if(APPLE)
      set(_pgo_wheel_runtime_path "@loader_path/../..")
    elseif(UNIX)
      set(_pgo_wheel_runtime_path "$ORIGIN/../..")
    endif()

    if(DEFINED _pgo_wheel_runtime_path)
      set_target_properties(${target} PROPERTIES
        BUILD_RPATH "${_pgo_wheel_runtime_path}"
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH "${_pgo_wheel_runtime_path}")
    endif()
    return()
  endif()

  if(TARGET pgo_windows_runtime)
    add_dependencies(${target} pgo_windows_runtime)
  endif()
  if(TARGET pgo_windows_runtime AND PGO_ATTACH_COPY_RUNTIME)
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${PGO_WINDOWS_RUNTIME_DLLS}
        "$<TARGET_FILE_DIR:${target}>"
      COMMENT "Staging Windows runtime DLLs beside $<TARGET_FILE_NAME:${target}>"
      COMMAND_EXPAND_LISTS
      VERBATIM)
  endif()
endfunction()
