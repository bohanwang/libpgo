cmake_minimum_required(VERSION 3.28.0 FATAL_ERROR)

include(CMakeModules/projectVersion.cmake)
pgo_read_project_version(
  "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml"
  PGO_PROJECT_VERSION)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile_commands.json for clangd" FORCE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(PGO_USE_MKL "Use MKL" ON)
option(PGO_ENABLE_PYTHON "Build python module" OFF)
option(PGO_ENABLE_FULL "Enable all functionalities" OFF)
option(PGO_ENABLE_CUDA "Enable CUDA" OFF)
option(PGO_BUILD_SUBPROJECTS "Include subprojects" OFF)
option(PGO_BUILD_TESTS "Build libpgo tests" OFF)

option(PGO_UI "Build pgo ui" OFF)
option(PGO_RENDER_USE_ASSIMP "Use Assimp" OFF)
option(PGO_RELEASE_MODE_DEBUG "Enable debug in release mode" OFF)
option(PGO_IGNORE_DEBUG_FLAG "Ignore debug flag" OFF)
option(PGO_CREATE_MISSING_FILE "Create file when it is missing" OFF)
option(PGO_HAS_ORIG_PARDISO "Use original pardiso" OFF)
option(PGO_OPT_USE_IPOPT "Use IPOPT" OFF)
option(PGO_OPT_USE_KNITRO "Use KNITRO" OFF)
option(PGO_MKL_LINK_DYNAMIC "MKL Link" ON)
option(PGO_ENABLE_ALEMBIC "Enable Alembic support" OFF)
option(PGO_ENABLE_GMSH "Enable GMSH support" OFF)
option(PGO_PORTABLE_BUILD "Disable host-specific code generation and release debug info" OFF)

if(APPLE AND PGO_ENABLE_PYTHON AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "26.0" CACHE STRING
    "Minimum macOS version for Python extension builds")
endif()

project(libpgo VERSION "${PGO_PROJECT_VERSION}" LANGUAGES CXX C)

if(NOT DEFINED PYPGO_VERSION_INFO)
  set(PYPGO_VERSION_INFO "${PROJECT_VERSION}")
endif()

if(APPLE)
  set(PGO_USE_MKL OFF CACHE BOOL "Use MKL" FORCE)
  message(STATUS "PGO_USE_MKL is forced OFF on macOS; using Accelerate")
endif()

# Keep the conventional CMake switch disabled for the entire configure so
# FetchContent dependencies cannot interpret it as permission to build their
# own test suites. PGO_BUILD_TESTS controls only libpgo's tests below.
set(BUILD_TESTING OFF CACHE BOOL "Third-party tests remain disabled" FORCE)

# Python packages such as mkl-devel install their CMake package files below the
# active Python prefix. Prefer the active virtual environment because `uv run`
# can launch CMake without placing the venv's Python first on PATH on Windows.
# This keeps uv/pip and direct CMake builds free of dependency-location flags.
if(DEFINED ENV{VIRTUAL_ENV} AND IS_DIRECTORY "$ENV{VIRTUAL_ENV}")
  set(_pgo_python_prefix "$ENV{VIRTUAL_ENV}")
elseif(PYTHON_EXECUTABLE)
  set(_pgo_python_executable "${PYTHON_EXECUTABLE}")
else()
  find_program(_pgo_python_executable NAMES python3 python)
endif()
if(NOT _pgo_python_prefix AND _pgo_python_executable)
  execute_process(
    COMMAND "${_pgo_python_executable}" -c "import sys; print(sys.prefix)"
    OUTPUT_VARIABLE _pgo_python_prefix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pgo_python_prefix_result
  )
  if(NOT _pgo_python_prefix_result EQUAL 0)
    unset(_pgo_python_prefix)
  endif()
endif()
if(_pgo_python_prefix AND IS_DIRECTORY "${_pgo_python_prefix}")
  list(PREPEND CMAKE_PREFIX_PATH "${_pgo_python_prefix}")
  if(WIN32 AND IS_DIRECTORY "${_pgo_python_prefix}/Library")
    list(PREPEND CMAKE_PREFIX_PATH "${_pgo_python_prefix}/Library")
  endif()
  set(PGO_PYTHON_DEPENDENCY_PREFIX "${_pgo_python_prefix}" CACHE INTERNAL
    "Python environment that provides native build dependencies" FORCE)
  message(STATUS "Using Python dependency prefix: ${_pgo_python_prefix}")
endif()
unset(_pgo_python_executable)
unset(_pgo_python_prefix)
unset(_pgo_python_prefix_result)

if(APPLE)
  set(BLA_VENDOR Apple CACHE STRING "BLAS/LAPACK provider" FORCE)
  set(BLA_SIZEOF_INTEGER 4 CACHE STRING "LP64 BLAS/LAPACK integer size" FORCE)
endif()

set(PGO_CLANGD_BUILD_DIR "${CMAKE_SOURCE_DIR}/build")
set(PGO_CLANGD_COMPILE_COMMANDS "${PGO_CLANGD_BUILD_DIR}/compile_commands.json")

add_custom_target(copy_compile_commands_json
  COMMAND ${CMAKE_COMMAND} -E make_directory "${PGO_CLANGD_BUILD_DIR}"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_BINARY_DIR}/compile_commands.json"
          "${PGO_CLANGD_COMPILE_COMMANDS}"
  BYPRODUCTS "${PGO_CLANGD_COMPILE_COMMANDS}"
  COMMENT "Copying compile_commands.json to ${PGO_CLANGD_COMPILE_COMMANDS}"
  VERBATIM
)

if (PGO_ENABLE_CUDA)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(EXISTS "/usr/bin/g++")
      message(STATUS "Set CUDA host compiler: /usr/bin/g++")
      set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++)
    else()
    endif()
  endif()

  set(CMAKE_CUDA_ARCHITECTURES native)
  enable_language(CUDA)
endif()

# message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")
set(PGO_CMAKE_MODULES ${CMAKE_SOURCE_DIR}/CMakeModules)
list(APPEND CMAKE_MODULE_PATH "${PGO_CMAKE_MODULES}")
message(STATUS "CMake Module Path: ${CMAKE_MODULE_PATH}")

include(CompilationUtils)

if(APPLE)
  find_program(
    _pgo_homebrew_executable
    NAMES brew
    HINTS /opt/homebrew/bin /usr/local/bin
    NO_CACHE
  )
  if(NOT _pgo_homebrew_executable)
    message(FATAL_ERROR
      "Homebrew is required to locate OpenMP on macOS. "
      "Install Homebrew, then run: brew install libomp")
  endif()

  execute_process(
    COMMAND "${_pgo_homebrew_executable}" --prefix libomp
    OUTPUT_VARIABLE _pgo_homebrew_libomp_prefix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pgo_homebrew_libomp_result
    ERROR_QUIET
  )
  if(NOT _pgo_homebrew_libomp_result EQUAL 0
     OR NOT EXISTS "${_pgo_homebrew_libomp_prefix}/include/omp.h"
     OR NOT EXISTS "${_pgo_homebrew_libomp_prefix}/lib/libomp.dylib")
    message(FATAL_ERROR
      "Homebrew libomp is required on macOS. Run: brew install libomp")
  endif()

  # libomp is keg-only, so it is not linked into Homebrew's common include
  # and library directories. Prefer its exact paths over stale Conda or other
  # OpenMP cache entries from an earlier configure.
  list(PREPEND CMAKE_PREFIX_PATH "${_pgo_homebrew_libomp_prefix}")
  set(OpenMP_CXX_INCLUDE_DIR
    "${_pgo_homebrew_libomp_prefix}/include"
    CACHE PATH "AppleClang OpenMP include directory" FORCE)
  set(OpenMP_libomp_LIBRARY
    "${_pgo_homebrew_libomp_prefix}/lib/libomp.dylib"
    CACHE FILEPATH "AppleClang OpenMP runtime" FORCE)
  message(STATUS "Using Homebrew libomp: ${_pgo_homebrew_libomp_prefix}")

  unset(_pgo_homebrew_executable)
  unset(_pgo_homebrew_libomp_prefix)
  unset(_pgo_homebrew_libomp_result)
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
message(STATUS "OpenMP CXX found: ${OpenMP_CXX_FLAGS}")

# CPU feature probes execute on the build host, so they are valid only for
# explicitly host-optimized builds and must never affect portable wheels.
if(PGO_PORTABLE_BUILD)
  message(STATUS "Portable build: disabling host CPU feature detection")
else()
  include(Find_AVX)
  CHECK_FOR_AVX()
endif()

# basic debug flags
add_library(compilation_flag_for_debug INTERFACE)
add_flag(compilation_flag_for_debug MSVC Release INTERFACE /ZI)
add_flag(compilation_flag_for_debug MSVC Release INTERFACE /Oi-)
add_flag(compilation_flag_for_debug MSVC Release INTERFACE /Ob0)
add_def(compilation_flag_for_debug MSVC Debug INTERFACE _DEBUG)
add_def(compilation_flag_for_debug MSVC Release INTERFACE NDEBUG)

add_flag_poxis(compilation_flag_for_debug INTERFACE -ggdb3)
add_flag_poxis(compilation_flag_for_debug INTERFACE -O0)
add_flag_poxis(compilation_flag_for_debug INTERFACE -fsanitize=address)
add_flag_poxis(compilation_flag_for_debug INTERFACE -fsanitize=undefined)

add_flag(compilation_flag_for_debug GNU All INTERFACE -fsanitize=leak)
add_flag(compilation_flag_for_debug Clang All INTERFACE -fsanitize=leak)

add_library(compilation_flag INTERFACE)
add_library(cuda_compilation_flag INTERFACE)

if(TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(compilation_flag INTERFACE OpenMP::OpenMP_CXX)
  target_compile_definitions(compilation_flag INTERFACE USE_OPENMP)
endif()

message(STATUS "PGO compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
  message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

  if(CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "CMake build type is empty and will be set to Release")
    set(CMAKE_BUILD_TYPE Release)
  endif()

  if(CMAKE_GENERATOR STREQUAL "Ninja")
    add_flag_poxis(compilation_flag INTERFACE -fdiagnostics-color=always)
  endif()

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  message(STATUS "Windows platform...")
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")

  if(CMAKE_GENERATOR STREQUAL Ninja)
    message(STATUS "use ninja")
  else()
    set(CMAKE_CONFIGURATION_TYPES Debug Release)
  endif()
else()
  message(FATAL_ERROR "Unsupported compiler")
endif()

# basic compilation flags
add_flag_poxis(compilation_flag INTERFACE -Wall)
add_flag_poxis(compilation_flag INTERFACE -Wextra)
add_flag_poxis(compilation_flag INTERFACE -frounding-math)
add_flag_poxis(compilation_flag INTERFACE -fvisibility=hidden)
add_flag(compilation_flag GNU All INTERFACE -Wno-template-body)

if(NOT PGO_PORTABLE_BUILD)
  add_flag_poxis(compilation_flag INTERFACE -march=native)
  add_flag_poxis(compilation_flag INTERFACE -mtune=native)
  add_link_flag_poxis(compilation_flag -march=native)
  add_link_flag_poxis(compilation_flag -mtune=native)
endif()

# static gcc
add_flag(compilation_flag GNU All INTERFACE -static-libstdc++)
add_flag(compilation_flag GNU All INTERFACE -static-libgcc)
add_link_flag(compilation_flag GNU All -static-libstdc++)
add_link_flag(compilation_flag GNU All -static-libgcc)

# msvc flags
add_def(compilation_flag MSVC All INTERFACE __BASE_FILE__=__FILE__)
add_flag(compilation_flag MSVC All INTERFACE /bigobj)
add_flag(compilation_flag MSVC All INTERFACE /Zc:__cplusplus)
add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/bigobj")

if(NOT PGO_PORTABLE_BUILD)
  if(HAVE_AVX512_EXTENSIONS)
    add_flag(compilation_flag MSVC All INTERFACE /arch:AVX512)
    add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/arch:AVX512")
  elseif(HAVE_AVX2_EXTENSIONS)
    add_flag(compilation_flag MSVC All INTERFACE /arch:AVX2)
    add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/arch:AVX2")
  elseif(HAVE_AVX_EXTENSIONS)
    add_flag(compilation_flag MSVC All INTERFACE /arch:AVX)
    add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/arch:AVX")
  endif()
endif()

add_flag(compilation_flag MSVC All INTERFACE /MP)
add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/MP")

# set(cuda_compilation_flag "${cuda_compilation_flag};/MP")
add_def(compilation_flag MSVC All INTERFACE NOMINMAX)
add_def(compilation_flag MSVC All INTERFACE _USE_MATH_DEFINES)
add_def(compilation_flag MSVC All INTERFACE _CRT_SECURE_NO_WARNINGS)
add_def(compilation_flag MSVC All INTERFACE HAVE_STRUCT_TIMESPEC)
add_def(compilation_flag MSVC All INTERFACE GLOG_NO_ABBREVIATED_SEVERITIES)

if(PGO_ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug Mode..")

  add_flag_poxis(compilation_flag INTERFACE -O0)
  add_flag_poxis(compilation_flag INTERFACE -ggdb3)
  add_link_flag_poxis(compilation_flag -ggdb3)
else()
  message(STATUS "Release Mode..")

  add_flag_poxis(compilation_flag INTERFACE -O3)
  if(NOT PGO_PORTABLE_BUILD)
    add_flag_poxis(compilation_flag INTERFACE -ggdb3)
    add_link_flag_poxis(compilation_flag -ggdb3)
  endif()
endif()

# ================== third party ==========================
find_package(Threads REQUIRED)

if(PGO_UI)
  find_package(Vulkan REQUIRED)

  if(PGO_RENDER_USE_ASSIMP)
    find_package(assimp REQUIRED)
  endif()

  include(${PROJECT_SOURCE_DIR}/CMakeModules/third-party/glfw.cmake)
  include(${PROJECT_SOURCE_DIR}/CMakeModules/third-party/glm.cmake)
endif()

add_subdirectory(third-party)

include(CMakeModules/third-party/tbb.cmake)
include(CMakeModules/third-party/mkl.cmake)
include(CMakeModules/third-party/backward.cmake)

if(PGO_ENABLE_PYTHON OR PGO_ENABLE_FULL)
  include(CMakeModules/third-party/suitesparse.cmake)
endif()

include(CMakeModules/third-party/eigen.cmake)
include(CMakeModules/third-party/fmt.cmake)
include(CMakeModules/third-party/nlohmann_json.cmake)
include(CMakeModules/third-party/spdlog.cmake)
include(CMakeModules/third-party/autodiff.cmake)
include(CMakeModules/third-party/argparse.cmake)

if(PGO_ENABLE_PYTHON OR PGO_ENABLE_FULL)
  include(CMakeModules/third-party/ceres.cmake)
  include(CMakeModules/third-party/boost.cmake)
  include(CMakeModules/third-party/cgal.cmake)
endif()

if(PGO_ENABLE_FULL)
  include(CMakeModules/third-party/geogram.cmake)
endif()

if(PGO_ENABLE_FULL OR PGO_ENABLE_PYTHON)
  include(CMakeModules/third-party/libigl.cmake)
  
  if(PGO_ENABLE_ALEMBIC)
    include(CMakeModules/third-party/imath.cmake)
    include(CMakeModules/third-party/alembic.cmake)
  endif()
  
  if(PGO_ENABLE_GMSH)
    include(CMakeModules/third-party/gmsh.cmake)
  endif()
endif()

if(PGO_ENABLE_PYTHON)
  include(CMakeModules/third-party/pybind11.cmake)
endif()

include(CMakeModules/windowsRuntime.cmake)

if(TARGET MKL::MKL)
  add_compile_definitions(PGO_HAS_MKL)
endif()

if(TARGET CGAL::CGAL)
  add_compile_definitions(PGO_HAS_CGAL)
endif()

if(TARGET Ceres::ceres)
  add_compile_definitions(PGO_HAS_CERES)
endif()

if(PGO_OPT_USE_KNITRO)
  include(CMakeModules/third-party/knitro.cmake)

  if(TARGET Knitro::Knitro)
    add_compile_definitions(PGO_HAS_KNITRO)
  endif()
endif()

if(PGO_HAS_ORIG_PARDISO)
  include(CMakeModules/third-party/pardiso.cmake)

  if(TARGET Pardiso::Pardiso)
    add_compile_definitions(PGO_HAS_ORIG_PARDISO)
  endif()
endif()

# =================== external ====================
if(MSVC)
  find_program(POWERSHELL_PATH powershell REQUIRED)
  message(STATUS "Powershell: ${POWERSHELL_PATH}")

  get_filename_component(VCPKG_TOOLCHAIN_FILE_DIR "${CMAKE_TOOLCHAIN_FILE}" DIRECTORY)
  set(APPLOCAL_FILE_PATH "${VCPKG_TOOLCHAIN_FILE_DIR}/msbuild/applocal.ps1")
  message(STATUS "applocal.ps1: ${APPLOCAL_FILE_PATH}")

  if(EXISTS ${APPLOCAL_FILE_PATH})
    set(PGO_HAS_APPLOCAL ON)
    message(STATUS "Found.")
  else()
    set(PGO_HAS_APPLOCAL OFF)
    message(STATUS "Not using vcpkg.")
  endif()
endif()

set(PGO_BINARY_ROOT ${CMAKE_CURRENT_BINARY_DIR})

function(add_libpgo_lib tgt sources headers)
  if(PGO_CREATE_MISSING_FILE)
    foreach(src ${sources})
      if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${src})
        message(STATUS "Create file ${CMAKE_CURRENT_SOURCE_DIR}/${src}")
        file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/${src} "")
      endif()
    endforeach()

    foreach(h ${headers})
      if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${h})
        message(STATUS "Create file ${CMAKE_CURRENT_SOURCE_DIR}/${h}")
        file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/${h} "")
      endif()
    endforeach()
  endif()

  add_library(${tgt} STATIC ${sources} ${headers})

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set_target_properties(${tgt} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PGO_BINARY_ROOT}/lib)
    set_target_properties(${tgt} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PGO_BINARY_ROOT}/lib)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(${tgt} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PGO_BINARY_ROOT}/lib/Debug)
    set_target_properties(${tgt} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PGO_BINARY_ROOT}/lib/Release)
  else()
    message(FATAL_ERROR "Unsupported compiler")
  endif()

  target_include_directories(${tgt} PUBLIC ./)
  target_link_libraries(${tgt} PRIVATE compilation_flag)

  if(PGO_RELEASE_MODE_DEBUG)
    target_link_libraries(${tgt} PRIVATE compilation_flag_for_debug)
  endif()

  set_property(TARGET ${tgt} PROPERTY FOLDER libraries)

  # set(PGO_GLOBAL_LIBRARY_TARGETS ${PGO_GLOBAL_LIBRARY_TARGETS} ${tgt} CACHE INTERNAL "global library targets")
endfunction()

function(add_libpgo_tools tgt sources deps)
  if(PGO_CREATE_MISSING_FILE)
    foreach(src ${sources})
      if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${src})
        message(STATUS "Create file ${CMAKE_CURRENT_SOURCE_DIR}/${src}")
        file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/${src} "")
      endif()
    endforeach()
  endif()

  add_executable(${tgt} ${sources})

  if(TARGET pgo_windows_runtime)
    add_dependencies(${tgt} pgo_windows_runtime)
  endif()

  foreach(d ${deps})
    target_link_libraries(${tgt} PRIVATE ${d})
  endforeach()

  fix_mkl_tbb_debug(${tgt})
  fix_glut32(${tgt})

  target_link_libraries(${tgt} PRIVATE compilation_flag)

  if(PGO_RELEASE_MODE_DEBUG)
    target_link_libraries(${tgt} PRIVATE compilation_flag_for_debug)
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set_target_properties(${tgt} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PGO_BINARY_ROOT}/bin)
    set_target_properties(${tgt} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PGO_BINARY_ROOT}/bin)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(${tgt} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PGO_BINARY_ROOT}/bin/Debug)
    set_target_properties(${tgt} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PGO_BINARY_ROOT}/bin/Release)
  endif()

  set_property(TARGET ${tgt} PROPERTY FOLDER "tools")

  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC" AND PGO_HAS_APPLOCAL)
    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${POWERSHELL_PATH} -noprofile -executionpolicy Bypass -file ${APPLOCAL_FILE_PATH}
      -targetBinary $<TARGET_FILE:${tgt}>
      -installedDir "${EXT_LIB_ROOT}/bin"
      -OutVariable out
    )

    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${POWERSHELL_PATH} -noprofile -executionpolicy Bypass -file ${APPLOCAL_FILE_PATH}
      -targetBinary $<TARGET_FILE:${tgt}>
      -installedDir "${EXT_LIB_ROOT}/bin/$<CONFIG>"
      -OutVariable out
    )
  endif()

  # if(PGO_COMPILER STREQUAL "GNU")
  # install(TARGETS ${tgt} RUNTIME DESTINATION bin)
  # elseif(PGO_COMPILER STREQUAL "MSVC")
  # install(TARGETS ${tgt} RUNTIME DESTINATION bin/$<CONFIG>)
  # endif()
endfunction()

add_subdirectory(src)

if(PGO_BUILD_TESTS)
  enable_testing()
  message(STATUS "libpgo testing enabled.")
  add_subdirectory(tests)
endif()

if(PGO_BUILD_SUBPROJECTS)
  add_subdirectory(projects)
endif()

if(PGO_ENABLE_PYTHON)
  include(CMakeModules/python_list.cmake)
endif()

# generate config file to be used in common build tree
configure_package_config_file(
  CMakeModules/pgoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/pgoConfig.cmake
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/"
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/pgoConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/pgoConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/pgoConfigVersion.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/"
)
