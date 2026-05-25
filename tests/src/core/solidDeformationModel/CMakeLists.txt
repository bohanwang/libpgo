add_executable(simulationMesh_gtest simulationMesh_gtest.cpp)
target_link_libraries(simulationMesh_gtest PRIVATE GTest::gtest_main solidDeformationModel volumetricMesh)
target_compile_definitions(simulationMesh_gtest PRIVATE
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/cubic/box/box.veg")
set_property(TARGET simulationMesh_gtest PROPERTY FOLDER "tests/gtest")
pgo_gtest_discover_tests(simulationMesh_gtest)

add_executable(deformationModelAssembler_gtest deformationModelAssembler_gtest.cpp)
target_link_libraries(deformationModelAssembler_gtest PRIVATE GTest::gtest_main solidDeformationModel volumetricMesh mesh)
target_compile_definitions(deformationModelAssembler_gtest PRIVATE
  LIBPGO_TEST_TORUS_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/tet/torus.veg"
  LIBPGO_TEST_SHELL_OBJ="${CMAKE_SOURCE_DIR}/tests/fixtures/shell/shell.obj"
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/cubic/box/box.veg")
set_property(TARGET deformationModelAssembler_gtest PROPERTY FOLDER "tests/gtest")
pgo_gtest_discover_tests(deformationModelAssembler_gtest)

add_executable(deformationModelFactory_gtest deformationModelFactory_gtest.cpp)
target_link_libraries(deformationModelFactory_gtest PRIVATE GTest::gtest_main solidDeformationModel volumetricMesh mesh)
target_compile_definitions(deformationModelFactory_gtest PRIVATE
  LIBPGO_TEST_TORUS_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/tet/torus.veg")
set_property(TARGET deformationModelFactory_gtest PROPERTY FOLDER "tests/gtest")
pgo_gtest_discover_tests(deformationModelFactory_gtest)

add_executable(cubicMeshDeformationModel_gtest cubicMeshDeformationModel_gtest.cpp)
target_link_libraries(cubicMeshDeformationModel_gtest PRIVATE GTest::gtest_main solidDeformationModel)
set_property(TARGET cubicMeshDeformationModel_gtest PROPERTY FOLDER "tests/gtest")
pgo_gtest_discover_tests(cubicMeshDeformationModel_gtest)

add_executable(deformationModelEnergyMaxStep_gtest deformationModelEnergyMaxStep_gtest.cpp)
target_link_libraries(deformationModelEnergyMaxStep_gtest PRIVATE GTest::gtest_main solidDeformationModel simulation mesh)
target_compile_definitions(deformationModelEnergyMaxStep_gtest PRIVATE
  LIBPGO_TEST_SHELL_OBJ="${CMAKE_SOURCE_DIR}/tests/fixtures/shell/shell.obj")
set_property(TARGET deformationModelEnergyMaxStep_gtest PROPERTY FOLDER "tests/gtest")
pgo_gtest_discover_tests(deformationModelEnergyMaxStep_gtest)
