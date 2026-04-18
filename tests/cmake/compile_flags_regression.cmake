cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED PROJECT_SOURCE_DIR)
  message(FATAL_ERROR "PROJECT_SOURCE_DIR is required")
endif()

if(NOT DEFINED TEST_ROOT)
  message(FATAL_ERROR "TEST_ROOT is required")
endif()

if(NOT DEFINED TEST_GENERATOR OR TEST_GENERATOR STREQUAL "")
  message(FATAL_ERROR "TEST_GENERATOR is required")
endif()

if(NOT DEFINED TEST_C_COMPILER OR TEST_C_COMPILER STREQUAL "")
  message(FATAL_ERROR "TEST_C_COMPILER is required")
endif()

if(NOT DEFINED TEST_CXX_COMPILER OR TEST_CXX_COMPILER STREQUAL "")
  message(FATAL_ERROR "TEST_CXX_COMPILER is required")
endif()

set(test_home "${TEST_ROOT}/home")
set(debug_build "${TEST_ROOT}/base_no_mkl_debug")
set(release_build "${TEST_ROOT}/base_no_mkl")
set(release_flags_make "${release_build}/src/core/contact/CMakeFiles/contact.dir/flags.make")
set(debug_eigen_export "${debug_build}/_deps/eigen3-build/Eigen3Targets.cmake")

file(REMOVE_RECURSE "${TEST_ROOT}")
file(MAKE_DIRECTORY "${test_home}")

function(run_nested_configure build_dir build_type)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env HOME=${test_home}
      "${CMAKE_COMMAND}"
      -S "${PROJECT_SOURCE_DIR}"
      -B "${build_dir}"
      -G "${TEST_GENERATOR}"
      -DCMAKE_C_COMPILER=${TEST_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${TEST_CXX_COMPILER}
      -DCMAKE_BUILD_TYPE=${build_type}
      -DBUILD_TESTING=OFF
      -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=FALSE
      -DCMAKE_FIND_USE_PACKAGE_REGISTRY=TRUE
      -DCMAKE_EXPORT_PACKAGE_REGISTRY=TRUE
      -DPGO_USE_MKL=OFF
      -DPGO_ENABLE_FULL=ON
      -DPGO_BUILD_SUBPROJECTS=ON
      -DPGO_ENABLE_ALEMBIC=ON
      -DPGO_ENABLE_GMSH=ON
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE output
  )

  if(NOT result EQUAL 0)
    message(FATAL_ERROR "Configure failed for ${build_type}:\n${output}")
  endif()
endfunction()

run_nested_configure("${debug_build}" "Debug")

if(NOT EXISTS "${debug_eigen_export}")
  message(FATAL_ERROR "Expected debug Eigen export at ${debug_eigen_export}")
endif()

file(READ "${debug_eigen_export}" debug_eigen_export_contents)
if(debug_eigen_export_contents MATCHES [[-O0]])
  message(FATAL_ERROR "Debug Eigen export leaked project compile options:\n${debug_eigen_export}")
endif()

run_nested_configure("${release_build}" "Release")

if(NOT EXISTS "${release_flags_make}")
  message(FATAL_ERROR "Expected release flags.make at ${release_flags_make}")
endif()

file(READ "${release_flags_make}" release_flags)

if(release_flags MATCHES [[(^|\n)CXX_FLAGS[^=\n]* = .* -O0([ \n]|$)]])
  message(FATAL_ERROR "Release flags unexpectedly contain -O0:\n${release_flags}")
endif()

if(NOT release_flags MATCHES [[-O3]])
  message(FATAL_ERROR "Release flags are missing -O3:\n${release_flags}")
endif()

if(NOT release_flags MATCHES [[-DNDEBUG]])
  message(FATAL_ERROR "Release flags are missing -DNDEBUG:\n${release_flags}")
endif()
