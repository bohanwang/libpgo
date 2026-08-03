add_executable(cubicMesher_gtest cubicMesher_gtest.cpp)
target_link_libraries(cubicMesher_gtest PRIVATE GTest::gtest_main cubicMesherCore)
target_include_directories(cubicMesher_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/cubicMesher)
add_dependencies(cubicMesher_gtest cubicMesher)
target_compile_definitions(cubicMesher_gtest PRIVATE
  PGO_TEST_CUBIC_MESHER_BIN="$<TARGET_FILE:cubicMesher>")
set_property(TARGET cubicMesher_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cubicMesher_gtest)

if(TARGET convertAnimation)
  add_executable(convertAnimation_gtest convertAnimation_gtest.cpp)
  target_link_libraries(convertAnimation_gtest PRIVATE GTest::gtest_main)
  add_dependencies(convertAnimation_gtest convertAnimation)
  target_compile_definitions(convertAnimation_gtest PRIVATE
    PGO_TEST_CONVERT_ANIMATION_BIN="$<TARGET_FILE:convertAnimation>")
  set_property(TARGET convertAnimation_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(convertAnimation_gtest)
endif()

add_executable(runSim_gtest
  runSim_gtest.cpp)
target_link_libraries(runSim_gtest PRIVATE
  GTest::gtest_main
  simulationRunner)
target_include_directories(runSim_gtest PRIVATE
  ${CMAKE_SOURCE_DIR}/src/simulationRunner
  ${CMAKE_SOURCE_DIR}/src/tools/runSim
  ${CMAKE_SOURCE_DIR}/src/core/solidDeformationModel)
target_compile_definitions(runSim_gtest PRIVATE
  LIBPGO_TEST_TET_BOX_VEG="${LIBPGO_TEST_TET_BOX_VEG}"
  LIBPGO_TEST_TET_BOX_OBJ="${LIBPGO_TEST_TET_BOX_OBJ}"
  LIBPGO_TEST_CUBIC_BOX_VEG="${LIBPGO_TEST_CUBIC_BOX_VEG}"
  LIBPGO_TEST_CUBIC_BOX_OBJ="${LIBPGO_TEST_CUBIC_BOX_OBJ}")
add_dependencies(runSim_gtest libpgo_test_cubic_assets)
set_property(TARGET runSim_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runSim_gtest)

add_executable(runIPCSim_gtest
  runIPCSim_gtest.cpp)
target_link_libraries(runIPCSim_gtest PRIVATE
  GTest::gtest_main
  simulationRunner)
target_include_directories(runIPCSim_gtest PRIVATE
  ${CMAKE_SOURCE_DIR}/src/simulationRunner
  ${CMAKE_SOURCE_DIR}/src/tools/runSim
  ${CMAKE_SOURCE_DIR}/src/core/basicIO
  ${CMAKE_SOURCE_DIR}/src/core/contact
  ${CMAKE_SOURCE_DIR}/src/core/constraintPotentialEnergies
  ${CMAKE_SOURCE_DIR}/src/core/eigenSupport
  ${CMAKE_SOURCE_DIR}/src/core/interpolationCoordinates
  ${CMAKE_SOURCE_DIR}/src/core/libiglInterface
  ${CMAKE_SOURCE_DIR}/src/core/mesh
  ${CMAKE_SOURCE_DIR}/src/core/pgoLogging
  ${CMAKE_SOURCE_DIR}/src/core/simulation
  ${CMAKE_SOURCE_DIR}/src/core/solidDeformationModel
  ${CMAKE_SOURCE_DIR}/src/core/volumetricMesh)
target_compile_definitions(runIPCSim_gtest PRIVATE
  LIBPGO_TEST_SHELL_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/examples/assets/shell"
  LIBPGO_TEST_IPC_TET_EXAMPLE_DIR="${LIBPGO_TEST_TET_BOX_DIR}"
  LIBPGO_TEST_IPC_CUBIC_EXAMPLE_DIR="${LIBPGO_TEST_CUBIC_BOX_DIR}"
  PGO_TEST_RUN_IPC_SIM_BIN="$<TARGET_FILE:runIPCSim>")
add_dependencies(runIPCSim_gtest runIPCSim libpgo_test_cubic_assets)
set_property(TARGET runIPCSim_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runIPCSim_gtest)
