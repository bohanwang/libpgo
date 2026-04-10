include_guard(GLOBAL)

# libpgo uses a single Conan-based dependency workflow for both native and
# Python builds. This bootstrap file is the required entry point for resolving
# C++ dependencies during configure.

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# The variables in _LIBPGO_TRY_COMPILE_VARS are used in CMake try_compile calls during configure, so they must be added to CMAKE_TRY_COMPILE_PLATFORM_VARIABLES to be visible in that context.
set(_LIBPGO_TRY_COMPILE_VARS
  CMAKE_BUILD_TYPE
  PGO_CONAN_OUTPUT_DIR
  PGO_CONAN_GENERATORS_DIR
  PGO_CONAN_CPPSTD
  PGO_CONAN_PROFILE_HOST
  PGO_CONAN_PROFILE_BUILD
  PGO_PROFILE_FULL
  PGO_FEATURE_PYTHON
  PGO_FEATURE_ANIMATION_IO
  PGO_FEATURE_GEOMETRY_STACK
  PGO_FEATURE_LIBIGL
  PGO_FEATURE_GEOGRAM
  PGO_FEATURE_GMSH
  PGO_FEATURE_MKL
  PGO_FEATURE_ARPACK
)

foreach(_LIBPGO_VAR IN LISTS _LIBPGO_TRY_COMPILE_VARS)
  if(NOT _LIBPGO_VAR IN_LIST CMAKE_TRY_COMPILE_PLATFORM_VARIABLES)
    list(APPEND CMAKE_TRY_COMPILE_PLATFORM_VARIABLES "${_LIBPGO_VAR}")
  endif()
endforeach()

# Convert a boolean CMake variable to "True"/"False" string for Conan command line
function(_libpgo_bool_to_conan out_var value)
  if(${value})
    set(${out_var} "True" PARENT_SCOPE)
  else()
    set(${out_var} "False" PARENT_SCOPE)
  endif()
endfunction()

# Finds the host Python interpreter
function(_libpgo_find_host_python out_var)
  if(DEFINED ENV{UV_PYTHON} AND EXISTS "$ENV{UV_PYTHON}")
    set(${out_var} "$ENV{UV_PYTHON}" PARENT_SCOPE)
    return()
  endif()

  find_program(_LIBPGO_HOST_PYTHON NAMES python3 python REQUIRED)
  set(${out_var} "${_LIBPGO_HOST_PYTHON}" PARENT_SCOPE)
endfunction()

# Computes the effective feature set based on the full profile and inter-feature dependencies, then updates the PGO_FEATURE_* variables accordingly. 
# This ensures that the Conan install command is always invoked with a consistent and complete set of features.
function(_libpgo_compute_effective_features)
  set(_LIBPGO_EFFECTIVE_PROFILE_FULL ${PGO_PROFILE_FULL})
  set(_LIBPGO_EFFECTIVE_PYTHON ${PGO_FEATURE_PYTHON})
  set(_LIBPGO_EFFECTIVE_ANIMATION_IO ${PGO_FEATURE_ANIMATION_IO})
  set(_LIBPGO_EFFECTIVE_GEOMETRY_STACK ${PGO_FEATURE_GEOMETRY_STACK})
  set(_LIBPGO_EFFECTIVE_LIBIGL ${PGO_FEATURE_LIBIGL})
  set(_LIBPGO_EFFECTIVE_GEOGRAM ${PGO_FEATURE_GEOGRAM})
  set(_LIBPGO_EFFECTIVE_GMSH ${PGO_FEATURE_GMSH})
  set(_LIBPGO_EFFECTIVE_MKL ${PGO_FEATURE_MKL})
  set(_LIBPGO_EFFECTIVE_ARPACK ${PGO_FEATURE_ARPACK})

  if(_LIBPGO_EFFECTIVE_PROFILE_FULL)
    set(_LIBPGO_EFFECTIVE_PYTHON ON)
    set(_LIBPGO_EFFECTIVE_ANIMATION_IO ON)
    set(_LIBPGO_EFFECTIVE_GEOMETRY_STACK ON)
    set(_LIBPGO_EFFECTIVE_LIBIGL ON)
    set(_LIBPGO_EFFECTIVE_GEOGRAM ON)
    set(_LIBPGO_EFFECTIVE_GMSH ON)
  endif()

  if(_LIBPGO_EFFECTIVE_PYTHON OR _LIBPGO_EFFECTIVE_LIBIGL OR _LIBPGO_EFFECTIVE_GEOGRAM)
    set(_LIBPGO_EFFECTIVE_GEOMETRY_STACK ON)
  endif()

  set(PGO_FEATURE_PYTHON "${_LIBPGO_EFFECTIVE_PYTHON}" CACHE BOOL "Build Python bindings" FORCE)
  set(PGO_FEATURE_ANIMATION_IO "${_LIBPGO_EFFECTIVE_ANIMATION_IO}" CACHE BOOL "Enable Alembic/Imath animation I/O" FORCE)
  set(PGO_FEATURE_GEOMETRY_STACK "${_LIBPGO_EFFECTIVE_GEOMETRY_STACK}" CACHE BOOL "Enable Boost/CGAL/Ceres/NLopt/SuiteSparse stack" FORCE)
  set(PGO_FEATURE_LIBIGL "${_LIBPGO_EFFECTIVE_LIBIGL}" CACHE BOOL "Enable libigl integration" FORCE)
  set(PGO_FEATURE_GEOGRAM "${_LIBPGO_EFFECTIVE_GEOGRAM}" CACHE BOOL "Enable geogram integration" FORCE)
  set(PGO_FEATURE_GMSH "${_LIBPGO_EFFECTIVE_GMSH}" CACHE BOOL "Enable Gmsh support" FORCE)
  set(PGO_FEATURE_MKL "${_LIBPGO_EFFECTIVE_MKL}" CACHE BOOL "Enable MKL support" FORCE)
  set(PGO_FEATURE_ARPACK "${_LIBPGO_EFFECTIVE_ARPACK}" CACHE BOOL "Enable ARPACK support" FORCE)
endfunction()

get_filename_component(_LIBPGO_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
_libpgo_compute_effective_features()

if(NOT DEFINED PGO_CONAN_CPPSTD OR PGO_CONAN_CPPSTD STREQUAL "")
  set(PGO_CONAN_CPPSTD "20" CACHE STRING "Required C++ standard passed to Conan settings.compiler.cppstd" FORCE)
endif()

set(_LIBPGO_CONAN_PROFILE_HOST "")
if(DEFINED PGO_CONAN_PROFILE_HOST AND NOT PGO_CONAN_PROFILE_HOST STREQUAL "")
  if(IS_ABSOLUTE "${PGO_CONAN_PROFILE_HOST}")
    set(_LIBPGO_CONAN_PROFILE_HOST "${PGO_CONAN_PROFILE_HOST}")
  else()
    get_filename_component(_LIBPGO_CONAN_PROFILE_HOST "${_LIBPGO_SOURCE_DIR}/${PGO_CONAN_PROFILE_HOST}" ABSOLUTE)
  endif()
endif()

set(_LIBPGO_CONAN_PROFILE_BUILD "")
if(DEFINED PGO_CONAN_PROFILE_BUILD AND NOT PGO_CONAN_PROFILE_BUILD STREQUAL "")
  if(IS_ABSOLUTE "${PGO_CONAN_PROFILE_BUILD}")
    set(_LIBPGO_CONAN_PROFILE_BUILD "${PGO_CONAN_PROFILE_BUILD}")
  else()
    get_filename_component(_LIBPGO_CONAN_PROFILE_BUILD "${_LIBPGO_SOURCE_DIR}/${PGO_CONAN_PROFILE_BUILD}" ABSOLUTE)
  endif()
endif()

if(DEFINED PGO_CONAN_GENERATORS_DIR AND NOT PGO_CONAN_GENERATORS_DIR STREQUAL "")
  get_filename_component(_LIBPGO_GENERATORS_DIR "${PGO_CONAN_GENERATORS_DIR}" ABSOLUTE)
  if(DEFINED PGO_CONAN_OUTPUT_DIR AND NOT PGO_CONAN_OUTPUT_DIR STREQUAL "")
    get_filename_component(_LIBPGO_CONAN_ROOT "${PGO_CONAN_OUTPUT_DIR}" ABSOLUTE)
  else()
    get_filename_component(_LIBPGO_CONAN_ROOT "${_LIBPGO_GENERATORS_DIR}" DIRECTORY)
    get_filename_component(_LIBPGO_CONAN_ROOT "${_LIBPGO_CONAN_ROOT}" DIRECTORY)
    get_filename_component(_LIBPGO_CONAN_ROOT "${_LIBPGO_CONAN_ROOT}" DIRECTORY)
  endif()
else()
  set(_LIBPGO_CONAN_ROOT "${CMAKE_BINARY_DIR}/.conan-build/${CMAKE_BUILD_TYPE}")
  set(_LIBPGO_GENERATORS_DIR "${_LIBPGO_CONAN_ROOT}/build/${CMAKE_BUILD_TYPE}/generators")
  set(PGO_CONAN_OUTPUT_DIR "${_LIBPGO_CONAN_ROOT}" CACHE PATH "Required libpgo Conan output directory" FORCE)
  set(PGO_CONAN_GENERATORS_DIR "${_LIBPGO_GENERATORS_DIR}" CACHE PATH "Required libpgo Conan generators directory" FORCE)
endif()

set(_LIBPGO_TOOLCHAIN_FILE "${_LIBPGO_GENERATORS_DIR}/conan_toolchain.cmake")
set(_LIBPGO_SIGNATURE_FILE "${_LIBPGO_CONAN_ROOT}/feature-signature.txt")

_libpgo_bool_to_conan(_LIBPGO_WITH_PYTHON PGO_FEATURE_PYTHON)
_libpgo_bool_to_conan(_LIBPGO_WITH_ANIMATION_IO PGO_FEATURE_ANIMATION_IO)
_libpgo_bool_to_conan(_LIBPGO_WITH_GEOMETRY_STACK PGO_FEATURE_GEOMETRY_STACK)
_libpgo_bool_to_conan(_LIBPGO_WITH_LIBIGL PGO_FEATURE_LIBIGL)
_libpgo_bool_to_conan(_LIBPGO_WITH_GEOGRAM PGO_FEATURE_GEOGRAM)
_libpgo_bool_to_conan(_LIBPGO_WITH_GMSH PGO_FEATURE_GMSH)
_libpgo_bool_to_conan(_LIBPGO_WITH_MKL PGO_FEATURE_MKL)
_libpgo_bool_to_conan(_LIBPGO_WITH_ARPACK PGO_FEATURE_ARPACK)
_libpgo_bool_to_conan(_LIBPGO_PROFILE_FULL PGO_PROFILE_FULL)

# Create a signature string based on the current configuration to determine if Conan install needs to be re-run. This avoids unnecessary Conan calls if the relevant configuration hasn't changed since the last install.
string(CONCAT _LIBPGO_SIGNATURE_CONTENT
  "build_type=${CMAKE_BUILD_TYPE}\n"
  "cppstd=${PGO_CONAN_CPPSTD}\n"
  "profile_host=${_LIBPGO_CONAN_PROFILE_HOST}\n"
  "profile_build=${_LIBPGO_CONAN_PROFILE_BUILD}\n"
  "profile_full=${_LIBPGO_PROFILE_FULL}\n"
  "with_python=${_LIBPGO_WITH_PYTHON}\n"
  "with_animation_io=${_LIBPGO_WITH_ANIMATION_IO}\n"
  "with_geometry_stack=${_LIBPGO_WITH_GEOMETRY_STACK}\n"
  "with_libigl=${_LIBPGO_WITH_LIBIGL}\n"
  "with_geogram=${_LIBPGO_WITH_GEOGRAM}\n"
  "with_gmsh=${_LIBPGO_WITH_GMSH}\n"
  "with_mkl=${_LIBPGO_WITH_MKL}\n"
  "with_arpack=${_LIBPGO_WITH_ARPACK}\n"
)

message(STATUS "libpgo Conan configuration signature:\n${_LIBPGO_SIGNATURE_CONTENT}")

set(_LIBPGO_SHOULD_INSTALL TRUE)
if(EXISTS "${_LIBPGO_TOOLCHAIN_FILE}" AND EXISTS "${_LIBPGO_SIGNATURE_FILE}")
  file(READ "${_LIBPGO_SIGNATURE_FILE}" _LIBPGO_EXISTING_SIGNATURE)
  if(_LIBPGO_EXISTING_SIGNATURE STREQUAL _LIBPGO_SIGNATURE_CONTENT)
    set(_LIBPGO_SHOULD_INSTALL FALSE)
    message(STATUS "libpgo Conan configuration matches existing signature, skipping Conan install.")
  else()
    message(STATUS "libpgo Conan configuration differs from existing signature, will run Conan install to update dependencies.")
  endif()
endif()

# If the signature file doesn't exist or doesn't match the current configuration, we need to run Conan install to ensure dependencies are up to date. Otherwise, we can skip directly to including the toolchain and generators.
if(_LIBPGO_SHOULD_INSTALL)
  find_program(_LIBPGO_CONAN_COMMAND NAMES conan REQUIRED)
  _libpgo_find_host_python(_LIBPGO_HOST_PYTHON)

  # Before running Conan install, we need to ensure that the local Conan recipes are exported so that the install command can find them. This is done by running a Python script provided in the repository. If this step fails, we cannot proceed with Conan install, so we treat it as a fatal error.
  execute_process(
    COMMAND "${_LIBPGO_HOST_PYTHON}" "${_LIBPGO_SOURCE_DIR}/conan/recipes/export_recipes.py"
    WORKING_DIRECTORY "${_LIBPGO_SOURCE_DIR}"
    RESULT_VARIABLE _LIBPGO_EXPORT_RESULT
    COMMAND_ECHO STDOUT
  )
  if(NOT _LIBPGO_EXPORT_RESULT EQUAL 0)
    message(FATAL_ERROR "libpgo requires its local Conan recipes to be exported before configure. See output above.")
  endif()

  file(MAKE_DIRECTORY "${_LIBPGO_CONAN_ROOT}")

  # Construct the Conan install command with all necessary arguments and options based on the current configuration. This includes settings for build type, C++ standard, feature options, and any specified profiles. The command is built as a list to handle arguments with spaces correctly.
  set(_LIBPGO_CONAN_INSTALL_COMMAND
    "${_LIBPGO_CONAN_COMMAND}" install "${_LIBPGO_SOURCE_DIR}"
    --output-folder "${_LIBPGO_CONAN_ROOT}"
    --build missing
    -s "build_type=${CMAKE_BUILD_TYPE}"
    -s "compiler.cppstd=${PGO_CONAN_CPPSTD}"
    -o "&:with_python=${_LIBPGO_WITH_PYTHON}"
    -o "&:with_animation_io=${_LIBPGO_WITH_ANIMATION_IO}"
    -o "&:with_geometry_stack=${_LIBPGO_WITH_GEOMETRY_STACK}"
    -o "&:with_libigl=${_LIBPGO_WITH_LIBIGL}"
    -o "&:with_geogram=${_LIBPGO_WITH_GEOGRAM}"
    -o "&:with_gmsh=${_LIBPGO_WITH_GMSH}"
    -o "&:with_mkl=${_LIBPGO_WITH_MKL}"
    -o "&:with_arpack=${_LIBPGO_WITH_ARPACK}"
    -o "&:profile_full=${_LIBPGO_PROFILE_FULL}"
  )

  if(NOT _LIBPGO_CONAN_PROFILE_HOST STREQUAL "")
    list(APPEND _LIBPGO_CONAN_INSTALL_COMMAND -pr:h "${_LIBPGO_CONAN_PROFILE_HOST}")
    message(STATUS "libpgo Conan host profile: ${_LIBPGO_CONAN_PROFILE_HOST}")
  endif()

  if(NOT _LIBPGO_CONAN_PROFILE_BUILD STREQUAL "")
    list(APPEND _LIBPGO_CONAN_INSTALL_COMMAND -pr:b "${_LIBPGO_CONAN_PROFILE_BUILD}")
    message(STATUS "libpgo Conan build profile: ${_LIBPGO_CONAN_PROFILE_BUILD}")
  endif()

  execute_process(
    COMMAND ${_LIBPGO_CONAN_INSTALL_COMMAND}
    WORKING_DIRECTORY "${_LIBPGO_SOURCE_DIR}"
    RESULT_VARIABLE _LIBPGO_INSTALL_RESULT
    COMMAND_ECHO STDOUT
  )
  if(NOT _LIBPGO_INSTALL_RESULT EQUAL 0)
    message(
      FATAL_ERROR
      "libpgo requires Conan for dependency resolution. Ensure the 'conan' command is installed and working, "
      "then review the output above for the failing package."
    )
  endif()

  file(WRITE "${_LIBPGO_SIGNATURE_FILE}" "${_LIBPGO_SIGNATURE_CONTENT}")
endif()

# After ensuring Conan dependencies are up to date, include the generated toolchain file and add the generators directory to CMake paths if they exist. This allows the rest of the CMake configuration to use the resolved dependencies and any custom CMake modules provided by the generators.
if(EXISTS "${_LIBPGO_TOOLCHAIN_FILE}")
  include("${_LIBPGO_TOOLCHAIN_FILE}")
  message(STATUS "Included Conan toolchain file: ${_LIBPGO_TOOLCHAIN_FILE}")
endif()

# If the generators directory exists, prepend it to CMAKE_PREFIX_PATH and CMAKE_MODULE_PATH so that any find_package calls or module includes can locate the Conan-generated files.
if(EXISTS "${_LIBPGO_GENERATORS_DIR}")
message(STATUS "Adding Conan generators directory to CMake paths: ${_LIBPGO_GENERATORS_DIR}")
  list(PREPEND CMAKE_PREFIX_PATH "${_LIBPGO_GENERATORS_DIR}")
  list(PREPEND CMAKE_MODULE_PATH "${_LIBPGO_GENERATORS_DIR}")
endif()
