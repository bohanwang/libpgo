include(FetchContent)
include(GoogleTest)

set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest installation" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock build" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "Use shared CRT for GoogleTest on Windows" FORCE)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(googletest)

add_executable(pgo_c_gtest pgo_c_gtest.cpp)
target_include_directories(pgo_c_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/c)
target_link_libraries(pgo_c_gtest PRIVATE GTest::gtest_main pgo_c)
target_compile_definitions(pgo_c_gtest PRIVATE LIBPGO_TEST_TORUS_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/tet/torus.veg")
set_property(TARGET pgo_c_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(pgo_c_gtest)

add_executable(mesh_gtest mesh_gtest.cpp)
target_include_directories(mesh_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/core/mesh)
target_link_libraries(mesh_gtest PRIVATE GTest::gtest_main mesh)
set_property(TARGET mesh_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(mesh_gtest)

# This regression checks POSIX/GCC-Clang compile flags (-O0/-O3/-ggdb3/-DNDEBUG)
# and has no MSVC equivalent, so only register it for those toolchains.
if(NOT MSVC)
  add_test(
    NAME compile_flags_regression
    COMMAND ${CMAKE_COMMAND}
      -DPROJECT_SOURCE_DIR=${CMAKE_SOURCE_DIR}
      -DTEST_ROOT=${CMAKE_BINARY_DIR}/compile_flags_regression
      -DTEST_GENERATOR=${CMAKE_GENERATOR}
      -DTEST_C_COMPILER=${CMAKE_C_COMPILER}
      -DTEST_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compile_flags_regression.cmake
  )
endif()

add_subdirectory(src)
