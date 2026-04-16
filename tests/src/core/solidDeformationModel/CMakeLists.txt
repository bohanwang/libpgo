add_executable(simulationMesh_gtest simulationMesh_gtest.cpp)
target_link_libraries(simulationMesh_gtest PRIVATE GTest::gtest_main solidDeformationModel volumetricMesh)
target_compile_definitions(simulationMesh_gtest PRIVATE
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/examples/cubic/box/box.veg")
set_property(TARGET simulationMesh_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(simulationMesh_gtest)

add_executable(deformationModelAssembler_gtest deformationModelAssembler_gtest.cpp)
target_link_libraries(deformationModelAssembler_gtest PRIVATE GTest::gtest_main solidDeformationModel volumetricMesh mesh)
target_compile_definitions(deformationModelAssembler_gtest PRIVATE
  LIBPGO_TEST_TORUS_VEG="${CMAKE_SOURCE_DIR}/examples/torus.veg"
  LIBPGO_TEST_SHELL_OBJ="${CMAKE_SOURCE_DIR}/examples/shell/shell.obj")
set_property(TARGET deformationModelAssembler_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(deformationModelAssembler_gtest)
