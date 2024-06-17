cmake_minimum_required(VERSION 3.28.0 FATAL_ERROR)

project(libpgo LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

OPTION(PGO_RELEASE_MODE_DEBUG "Enable debug in release mode" OFF)
OPTION(PGO_IGNORE_DEBUG_FLAG "Ignore debug flag" OFF)
OPTION(PGO_CREATE_MISSING_FILE "Create file when it is missing" OFF)
OPTION(PGO_USE_MKL "Use MKL" OFF)
OPTION(PGO_ENABLE_FULL "Enable all functionalities" OFF)
OPTION(PGO_OPT_USE_IPOPT "Use IPOPT" OFF)
OPTION(PGO_OPT_USE_KNITRO "Use KNITRO" OFF)
option(PGO_ENABLE_PYTHON "Build python module" OFF)
option(PGO_MKL_LINK_DYNAMIC "MKL Link" ON)

set(PGO_CMAKE_MODULES ${CMAKE_SOURCE_DIR}/CMakeModules)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PGO_CMAKE_MODULES})
message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")

include(CompilationUtils)

find_package(OpenMP)

if(OPENMP_FOUND)
  set(OPENMP_FLAG "${OpenMP_CXX_FLAGS}")
  add_compile_definitions(USE_OPENMP)
endif()

include(Find_AVX)
CHECK_FOR_AVX()

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

include(CMakeModules/third-party/tbb.cmake)
include(CMakeModules/third-party/mkl.cmake)

if(PGO_ENABLE_FULL)
  include(CMakeModules/third-party/suitesparse.cmake)
endif()

include(CMakeModules/third-party/eigen.cmake)
include(CMakeModules/third-party/fmt.cmake)
include(CMakeModules/third-party/nlohmann_json.cmake)
include(CMakeModules/third-party/spdlog.cmake)
include(CMakeModules/third-party/autodiff.cmake)

if(PGO_OPT_USE_KNITRO)
  include(CMakeModules/third-party/knitro.cmake)

  if(TARGET Knitro::Knitro)
    add_compile_definitions(PGO_HAS_KNITRO)
  endif()
endif()

if(PGO_ENABLE_FULL)
  include(CMakeModules/third-party/ceres.cmake)
  include(CMakeModules/third-party/boost.cmake)
  include(CMakeModules/third-party/cgal.cmake)
  # include(CMakeModules/third-party/arpackng.cmake)
endif()

if (PGO_ENABLE_PYTHON)
  include(CMakeModules/third-party/pybind11.cmake)
endif()

if(TARGET MKL::MKL)
  add_compile_definitions(PGO_HAS_MKL)
endif()

if(TARGET CGAL::CGAL)
  add_compile_definitions(PGO_HAS_CGAL)
endif()

if(TARGET Ceres::ceres)
  add_compile_definitions(PGO_HAS_CERES)
else()
endif()

# =================== external ====================
add_subdirectory(third-party)

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

# generate config file to be used in common build tree
configure_package_config_file(
  CMakeModules/pgoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/pgoConfig.cmake
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/"
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/pgoConfig.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/"
)
