cmake_minimum_required(VERSION 3.28 FATAL_ERROR)

project(pythonDependenciesSmoke NONE)

foreach(_pgo_required
    LIBPGO_SOURCE_DIR
    EXPECTED_PYTHON_PREFIX)
  if(NOT DEFINED ${_pgo_required} OR "${${_pgo_required}}" STREQUAL "")
    message(FATAL_ERROR "${_pgo_required} is required")
  endif()
endforeach()

set(ENV{VIRTUAL_ENV} "${EXPECTED_PYTHON_PREFIX}")
set(PYTHON_EXECUTABLE "/stale/cmake/cache/python"
  CACHE FILEPATH "Synthetic stale Python" FORCE)
set(PGO_ENABLE_PYTHON ON)

include("${LIBPGO_SOURCE_DIR}/CMakeModules/pythonDependencies.cmake")
pgo_configure_python_dependencies()

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
    "import pathlib, sys; print(pathlib.Path(sys.prefix).resolve())"
  OUTPUT_VARIABLE _pgo_actual_prefix
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _pgo_python_result)

if(NOT _pgo_python_result EQUAL 0)
  message(FATAL_ERROR "Resolved Python executable cannot run")
endif()
if(NOT "${_pgo_actual_prefix}" STREQUAL "${EXPECTED_PYTHON_PREFIX}")
  message(FATAL_ERROR
    "Python ABI prefix ${_pgo_actual_prefix} does not match "
    "dependency prefix ${EXPECTED_PYTHON_PREFIX}")
endif()
if(NOT "${PGO_PYTHON_DEPENDENCY_PREFIX}" STREQUAL
   "${EXPECTED_PYTHON_PREFIX}")
  message(FATAL_ERROR
    "Unexpected dependency prefix: ${PGO_PYTHON_DEPENDENCY_PREFIX}")
endif()
