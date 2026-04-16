add_executable(cubicMesher_gtest cubicMesher_gtest.cpp)
target_link_libraries(cubicMesher_gtest PRIVATE GTest::gtest_main cubicMesherCore)
target_include_directories(cubicMesher_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/cubicMesher)
add_dependencies(cubicMesher_gtest cubicMesher)
target_compile_definitions(cubicMesher_gtest PRIVATE
  PGO_TEST_CUBIC_MESHER_BIN="$<TARGET_FILE:cubicMesher>")
set_property(TARGET cubicMesher_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cubicMesher_gtest)

add_executable(runSim_gtest
  runSim_gtest.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimVolumeMeshIO.cpp)
target_link_libraries(runSim_gtest PRIVATE
  GTest::gtest_main
  configFileJSON
  interpolationCoordinates
  volumetricMesh
  mesh)
target_include_directories(runSim_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/runSim)
target_compile_definitions(runSim_gtest PRIVATE
  LIBPGO_TEST_TET_BOX_VEG="${CMAKE_SOURCE_DIR}/examples/box/box.veg"
  LIBPGO_TEST_TET_BOX_OBJ="${CMAKE_SOURCE_DIR}/examples/box/box.obj"
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/examples/cubic/box/box.veg"
  LIBPGO_TEST_CUBIC_BOX_OBJ="${CMAKE_SOURCE_DIR}/examples/cubic/box/box.obj")
set_property(TARGET runSim_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runSim_gtest)
