include_guard(GLOBAL)

function(pgo_configure_python_dependencies)
  set(_pgo_python_executable)

  # An active virtual environment is the source of truth for both the Python
  # ABI and Python-distributed native dependencies. This intentionally replaces
  # a stale PYTHON_EXECUTABLE left in an existing CMake cache.
  if(DEFINED ENV{VIRTUAL_ENV} AND IS_DIRECTORY "$ENV{VIRTUAL_ENV}")
    if(WIN32)
      set(_pgo_venv_python "$ENV{VIRTUAL_ENV}/Scripts/python.exe")
    else()
      set(_pgo_venv_python "$ENV{VIRTUAL_ENV}/bin/python")
    endif()
    if(NOT EXISTS "${_pgo_venv_python}")
      message(FATAL_ERROR
        "VIRTUAL_ENV does not contain its expected Python executable: "
        "${_pgo_venv_python}")
    endif()
    set(_pgo_python_executable "${_pgo_venv_python}")
  elseif(PYTHON_EXECUTABLE)
    if(NOT EXISTS "${PYTHON_EXECUTABLE}")
      message(FATAL_ERROR
        "PYTHON_EXECUTABLE does not exist: ${PYTHON_EXECUTABLE}")
    endif()
    set(_pgo_python_executable "${PYTHON_EXECUTABLE}")
  else()
    find_program(_pgo_python_executable
      NAMES python3 python
      NO_CACHE)
  endif()

  # Never retain a dependency prefix from a previous Python environment.
  unset(PGO_PYTHON_DEPENDENCY_PREFIX CACHE)

  if(NOT _pgo_python_executable)
    if(PGO_ENABLE_PYTHON)
      message(FATAL_ERROR
        "PGO_ENABLE_PYTHON=ON requires a Python executable. Activate the "
        "project environment or set PYTHON_EXECUTABLE.")
    endif()
    return()
  endif()

  get_filename_component(
    _pgo_python_executable "${_pgo_python_executable}" ABSOLUTE)
  execute_process(
    COMMAND "${_pgo_python_executable}" -c
      "import pathlib, sys; print(pathlib.Path(sys.prefix).resolve())"
    OUTPUT_VARIABLE _pgo_python_prefix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pgo_python_prefix_result)
  if(NOT _pgo_python_prefix_result EQUAL 0
     OR NOT IS_DIRECTORY "${_pgo_python_prefix}")
    message(FATAL_ERROR
      "Failed to query a valid sys.prefix from ${_pgo_python_executable}")
  endif()

  # pybind11 2.x consumes the legacy PYTHON_EXECUTABLE cache entry. Force it
  # to the same interpreter used to derive the native dependency prefix.
  set(PYTHON_EXECUTABLE "${_pgo_python_executable}"
    CACHE FILEPATH "Python executable used by libpgo and pybind11" FORCE)
  set(Python_EXECUTABLE "${_pgo_python_executable}"
    CACHE FILEPATH "Python executable used by CMake FindPython" FORCE)
  set(Python3_EXECUTABLE "${_pgo_python_executable}"
    CACHE FILEPATH "Python executable used by CMake FindPython3" FORCE)
  set(PGO_PYTHON_DEPENDENCY_PREFIX "${_pgo_python_prefix}"
    CACHE INTERNAL
    "Python environment that provides native build dependencies" FORCE)

  set(_pgo_prefix_path ${CMAKE_PREFIX_PATH})
  list(PREPEND _pgo_prefix_path "${_pgo_python_prefix}")
  if(WIN32 AND IS_DIRECTORY "${_pgo_python_prefix}/Library")
    list(PREPEND _pgo_prefix_path "${_pgo_python_prefix}/Library")
  endif()
  list(REMOVE_DUPLICATES _pgo_prefix_path)
  set(CMAKE_PREFIX_PATH "${_pgo_prefix_path}" PARENT_SCOPE)

  message(STATUS "Using Python executable: ${_pgo_python_executable}")
  message(STATUS "Using Python dependency prefix: ${_pgo_python_prefix}")
endfunction()
