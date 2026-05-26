if(TARGET glm::glm)
  return()
endif()

message(STATUS "Loading glm...")

pgo_dep_option(GLM_BUILD_LIBRARY BOOL OFF "Build dynamic/static library")
pgo_dep_option(GLM_BUILD_TESTS BOOL OFF "Build the test programs")
pgo_dep_option(GLM_BUILD_INSTALL BOOL OFF "Generate the install target")
pgo_dep_option(GLM_ENABLE_CXX_20 BOOL ON "Enable C++ 20")
pgo_dep_option(GLM_ENABLE_LANG_EXTENSIONS BOOL OFF "Enable language extensions")
pgo_dep_option(GLM_ENABLE_FAST_MATH BOOL OFF "Enable fast math optimizations")
pgo_dep_option(GLM_ENABLE_SIMD_AVX BOOL ON "Enable AVX optimizations")
pgo_dep_option(GLM_ENABLE_SIMD_AVX2 BOOL OFF "Enable AVX2 optimizations")

include(FetchContent)
FetchContent_Declare(
  glm
  URL https://github.com/g-truc/glm/archive/refs/tags/1.0.1.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(glm)

message(STATUS "Done.")
