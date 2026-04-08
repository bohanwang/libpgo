cmake_minimum_required(VERSION 3.28.0 FATAL_ERROR)

project(libpgo LANGUAGES CXX C)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ================= options ==========================
option(PGO_PROFILE_FULL "Enable the full libpgo feature profile" OFF)
option(PGO_FEATURE_PYTHON "Build Python bindings" OFF)
option(PGO_FEATURE_ANIMATION_IO "Enable Alembic/Imath animation I/O" OFF)
option(PGO_FEATURE_GEOMETRY_STACK "Enable Boost/CGAL/Ceres/NLopt/SuiteSparse stack" OFF)
option(PGO_FEATURE_LIBIGL "Enable libigl integration" OFF)
option(PGO_FEATURE_GEOGRAM "Enable geogram integration" OFF)
option(PGO_FEATURE_GMSH "Enable Gmsh support" OFF)
option(PGO_FEATURE_MKL "Enable MKL support" OFF)
option(PGO_FEATURE_ARPACK "Enable ARPACK support" OFF)
option(PGO_BUILD_SUBPROJECTS "Include subprojects" OFF)
option(PGO_ENABLE_TESTS "Build libpgo tests" ON)
set(PGO_MSVC_ARCH "default" CACHE STRING "MSVC target ISA for Windows builds: default, AVX, AVX2, AVX512")
set_property(CACHE PGO_MSVC_ARCH PROPERTY STRINGS default AVX AVX2 AVX512)

option(PGO_RELEASE_MODE_DEBUG "Enable debug in release mode" OFF)
option(PGO_IGNORE_DEBUG_FLAG "Ignore debug flag" OFF)
option(PGO_CREATE_MISSING_FILE "Create file when it is missing" OFF)
option(PGO_OPT_USE_IPOPT "Use IPOPT" OFF)
option(PGO_OPT_USE_KNITRO "Use KNITRO" OFF)
option(PGO_MKL_LINK_DYNAMIC "MKL Link" ON)

# resolve effective feature options based on profile and dependencies

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

set(PGO_PROFILE_FULL ${_LIBPGO_EFFECTIVE_PROFILE_FULL} CACHE BOOL "Enable the full libpgo feature profile" FORCE)
set(PGO_FEATURE_PYTHON ${_LIBPGO_EFFECTIVE_PYTHON} CACHE BOOL "Build Python bindings" FORCE)
set(PGO_FEATURE_ANIMATION_IO ${_LIBPGO_EFFECTIVE_ANIMATION_IO} CACHE BOOL "Enable Alembic/Imath animation I/O" FORCE)
set(PGO_FEATURE_GEOMETRY_STACK ${_LIBPGO_EFFECTIVE_GEOMETRY_STACK} CACHE BOOL "Enable Boost/CGAL/Ceres/NLopt/SuiteSparse stack" FORCE)
set(PGO_FEATURE_LIBIGL ${_LIBPGO_EFFECTIVE_LIBIGL} CACHE BOOL "Enable libigl integration" FORCE)
set(PGO_FEATURE_GEOGRAM ${_LIBPGO_EFFECTIVE_GEOGRAM} CACHE BOOL "Enable geogram integration" FORCE)
set(PGO_FEATURE_GMSH ${_LIBPGO_EFFECTIVE_GMSH} CACHE BOOL "Enable Gmsh support" FORCE)
set(PGO_FEATURE_MKL ${_LIBPGO_EFFECTIVE_MKL} CACHE BOOL "Enable MKL support" FORCE)
set(PGO_FEATURE_ARPACK ${_LIBPGO_EFFECTIVE_ARPACK} CACHE BOOL "Enable ARPACK support" FORCE)

if(PGO_ENABLE_TESTS)
  enable_testing()
endif()

set(PGO_CMAKE_MODULES ${CMAKE_SOURCE_DIR}/CMakeModules)
list(APPEND CMAKE_MODULE_PATH "${PGO_CMAKE_MODULES}")
message(STATUS "CMake Module Path: ${CMAKE_MODULE_PATH}")

# copy compile command to build directory for tools like clang-tidy and clangd
include(CompilationUtils)
copy_compile_command()

find_package(OpenMP)

if(OPENMP_FOUND)
  set(OPENMP_FLAG "${OpenMP_CXX_FLAGS}")
  add_compile_definitions(USE_OPENMP)
endif()

# openmp may cause mkl pardiso error
include(Find_AVX)
if(MSVC)
  message(STATUS "Skipping automatic /arch selection on MSVC; configure PGO_MSVC_ARCH to opt into AVX targets")
else()
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

message(STATUS "PGO compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  if(EXISTS "/usr/bin/g++")
    message(STATUS "Set CUDA host compiler: /usr/bin/g++")
    set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++)
  else()
  endif()
endif()

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
add_flag_poxis(compilation_flag INTERFACE -march=native)
add_flag_poxis(compilation_flag INTERFACE -mtune=native)
add_flag(compilation_flag GNU All INTERFACE -Wno-template-body)

add_link_flag_poxis(compilation_flag -march=native)
add_link_flag_poxis(compilation_flag -mtune=native)

# static gcc
add_flag(compilation_flag GNU All INTERFACE -static-libstdc++)
add_flag(compilation_flag GNU All INTERFACE -static-libgcc)
add_link_flag(compilation_flag GNU All -static-libstdc++)
add_link_flag(compilation_flag GNU All -static-libgcc)

# if has openmp flags
if(OPENMP_FOUND)
  add_flag(compilation_flag GNU All INTERFACE ${OPENMP_FLAG})
  add_flag(compilation_flag Clang All INTERFACE ${OPENMP_FLAG})
  add_link_flag(compilation_flag GNU All ${OPENMP_FLAG})
  add_link_flag(compilation_flag Clang All ${OPENMP_FLAG})
  add_flag(compilation_flag MSVC All INTERFACE ${OPENMP_FLAG})
endif()

# msvc flags
add_def(compilation_flag MSVC All INTERFACE __BASE_FILE__=__FILE__)
add_flag(compilation_flag MSVC All INTERFACE /bigobj)
add_flag(compilation_flag MSVC All INTERFACE /Zc:__cplusplus)
add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/bigobj")

if(MSVC)
  string(TOUPPER "${PGO_MSVC_ARCH}" PGO_MSVC_ARCH_UPPER)

  if(PGO_MSVC_ARCH_UPPER STREQUAL "DEFAULT")
    message(STATUS "MSVC target ISA: default")
  elseif(PGO_MSVC_ARCH_UPPER STREQUAL "AVX")
    add_flag(compilation_flag MSVC All INTERFACE /arch:AVX)
    add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/arch:AVX")
    message(STATUS "MSVC target ISA: AVX")
  elseif(PGO_MSVC_ARCH_UPPER STREQUAL "AVX2")
    add_flag(compilation_flag MSVC All INTERFACE /arch:AVX2)
    add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/arch:AVX2")
    message(STATUS "MSVC target ISA: AVX2")
  elseif(PGO_MSVC_ARCH_UPPER STREQUAL "AVX512")
    add_flag(compilation_flag MSVC All INTERFACE /arch:AVX512)
    add_flag_cuda(cuda_compilation_flag MSVC All INTERFACE "-Xcompiler=/arch:AVX512")
    message(STATUS "MSVC target ISA: AVX512")
  else()
    message(FATAL_ERROR "Unsupported PGO_MSVC_ARCH='${PGO_MSVC_ARCH}'. Expected one of: default, AVX, AVX2, AVX512")
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

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Debug Mode..")

  add_flag_poxis(compilation_flag INTERFACE -O0)
  add_flag_poxis(compilation_flag INTERFACE -ggdb3)
  add_link_flag_poxis(compilation_flag -ggdb3)
else()
  message(STATUS "Release Mode..")

  add_flag_poxis(compilation_flag INTERFACE -O3)
endif()

# ================== third party ==========================
find_package(Threads REQUIRED)

if(DEFINED PGO_CONAN_GENERATORS_DIR AND EXISTS "${PGO_CONAN_GENERATORS_DIR}")
  list(PREPEND CMAKE_PREFIX_PATH "${PGO_CONAN_GENERATORS_DIR}")
  list(PREPEND CMAKE_MODULE_PATH "${PGO_CONAN_GENERATORS_DIR}")
endif()

# use conan to manage dependencies. The conan configuration will be bootstrapped if not found, which may take a while on the first run.
include(${CMAKE_SOURCE_DIR}/cmake/LibpgoConanDeps.cmake)

# Older libigl integration used igl_include() to add copyleft modules. The
# Conan package we consume now exposes the umbrella targets directly and may not
# define igl_include(), so only call it when that helper exists.
if(PGO_FEATURE_LIBIGL AND COMMAND igl_include)
  igl_include(copyleft core)
  igl_include(copyleft cgal)
endif()


# MKL uses system detection (not managed by Conan).
include(CMakeModules/third-party/mkl.cmake)

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

  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(_LIBPGO_RUNTIME_DLL_MANIFEST "${CMAKE_CURRENT_BINARY_DIR}/${tgt}-runtime-dlls-$<CONFIG>.cmake")
    file(GENERATE
      OUTPUT "${_LIBPGO_RUNTIME_DLL_MANIFEST}"
      CONTENT "set(libpgo_target_runtime_dlls \"$<TARGET_RUNTIME_DLLS:${tgt}>\")\n"
    )

    add_custom_command(TARGET ${tgt} POST_BUILD
      COMMAND ${CMAKE_COMMAND}
      -Dlibpgo_runtime_manifest="${_LIBPGO_RUNTIME_DLL_MANIFEST}"
      -Dlibpgo_runtime_output_dir="$<TARGET_FILE_DIR:${tgt}>"
      -P "${CMAKE_SOURCE_DIR}/cmake/CopyRuntimeDllsFromManifest.cmake"
      VERBATIM
    )
  endif()

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
endfunction()

add_subdirectory(src)

if(PGO_ENABLE_TESTS)
  add_subdirectory(tests)
  message(STATUS "Tests enabled, added tests subdirectory")
endif()

if(PGO_BUILD_SUBPROJECTS)
  add_subdirectory(projects)
  message(STATUS "Subprojects enabled, added projects subdirectory")
endif()

if(PGO_FEATURE_PYTHON)
  include(CMakeModules/python_list.cmake)
  message(STATUS "Python bindings enabled, added CMakeModules/python_list.cmake")
endif()

if(NOT SKBUILD)
  # generate config file to be used in common build tree
  message(STATUS "Generating libpgoConfig.cmake in build directory for use in common build tree")
  include(CMakePackageConfigHelpers)
  configure_package_config_file(
    CMakeModules/pgoConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/pgoConfig.cmake
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/"
  )
  message(STATUS "Installing libpgoConfig.cmake to ${CMAKE_INSTALL_LIBDIR}/cmake/pgo/")
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/pgoConfig.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/"
  )
else()
  message(STATUS "SKBUILD detected, skipping generation of libpgoConfig.cmake")
endif()
