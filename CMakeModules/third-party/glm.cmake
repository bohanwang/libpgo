if(TARGET glm::glm)
  return()
endif()

message(STATUS "Loading glm...")

set(GLM_BUILD_LIBRARY OFF CACHE BOOL "Build dynamic/static library" FORCE)
set(GLM_BUILD_TESTS OFF CACHE BOOL "Build the test programs" FORCE)
set(GLM_BUILD_INSTALL OFF CACHE BOOL "Generate the install target" FORCE)
set(GLM_ENABLE_CXX_20 ON CACHE BOOL "Enable C++ 20" FORCE)
set(GLM_ENABLE_LANG_EXTENSIONS OFF CACHE BOOL "Enable language extensions" FORCE)
set(GLM_ENABLE_FAST_MATH OFF CACHE BOOL "Enable fast math optimizations" FORCE)
set(GLM_ENABLE_SIMD_AVX ON CACHE BOOL "Enable AVX optimizations" FORCE)
set(GLM_ENABLE_SIMD_AVX2 OFF CACHE BOOL "Enable AVX2 optimizations" FORCE)

include(FetchContent)
FetchContent_Declare(
  glm
  URL https://github.com/g-truc/glm/archive/refs/tags/1.0.1.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES glm
)

FetchContent_MakeAvailable(glm)

message(STATUS "Done.")
