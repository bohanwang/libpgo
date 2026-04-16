add_executable(simulationMesh_gtest simulationMesh_gtest.cpp)
target_link_libraries(simulationMesh_gtest PRIVATE GTest::gtest_main solidDeformationModel volumetricMesh)
target_compile_definitions(simulationMesh_gtest PRIVATE
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/examples/cubic/box/box.veg")
set_property(TARGET simulationMesh_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(simulationMesh_gtest)
