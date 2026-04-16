add_executable(cubicMesher_gtest cubicMesher_gtest.cpp)
target_link_libraries(cubicMesher_gtest PRIVATE GTest::gtest_main cubicMesherCore)
target_include_directories(cubicMesher_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/cubicMesher)
add_dependencies(cubicMesher_gtest cubicMesher)
target_compile_definitions(cubicMesher_gtest PRIVATE
  PGO_TEST_CUBIC_MESHER_BIN="$<TARGET_FILE:cubicMesher>")
set_property(TARGET cubicMesher_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cubicMesher_gtest)
