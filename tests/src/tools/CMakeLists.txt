add_executable(cubicMesher_gtest cubicMesher_gtest.cpp)
target_link_libraries(cubicMesher_gtest PRIVATE GTest::gtest_main cubicMesherCore)
target_include_directories(cubicMesher_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/cubicMesher)
add_dependencies(cubicMesher_gtest cubicMesher)
target_compile_definitions(cubicMesher_gtest PRIVATE
  PGO_TEST_CUBIC_MESHER_BIN="$<TARGET_FILE:cubicMesher>")
set_property(TARGET cubicMesher_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cubicMesher_gtest)

add_executable(hexJacobianValidation_gtest hexJacobianValidation_gtest.cpp)
target_link_libraries(hexJacobianValidation_gtest PRIVATE GTest::gtest_main)
target_include_directories(hexJacobianValidation_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/runSim)
set_property(TARGET hexJacobianValidation_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(hexJacobianValidation_gtest)

add_executable(triangleAreaValidation_gtest triangleAreaValidation_gtest.cpp)
target_link_libraries(triangleAreaValidation_gtest PRIVATE GTest::gtest_main)
target_include_directories(triangleAreaValidation_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/runSim)
set_property(TARGET triangleAreaValidation_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(triangleAreaValidation_gtest)

add_executable(validateIPCTrajectory_gtest validateIPCTrajectory_gtest.cpp)
target_link_libraries(validateIPCTrajectory_gtest PRIVATE GTest::gtest_main simulationRunner)
add_dependencies(validateIPCTrajectory_gtest validateIPCTrajectory)
target_compile_definitions(validateIPCTrajectory_gtest PRIVATE
  PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN="$<TARGET_FILE:validateIPCTrajectory>"
  PGO_TEST_VALIDATOR_BOX_DIR="${CMAKE_SOURCE_DIR}/examples/assets/volume/box")
set_property(TARGET validateIPCTrajectory_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(validateIPCTrajectory_gtest)

if(TARGET convertAnimation)
  add_executable(convertAnimation_gtest convertAnimation_gtest.cpp)
  target_link_libraries(convertAnimation_gtest PRIVATE GTest::gtest_main)
  add_dependencies(convertAnimation_gtest convertAnimation)
  target_compile_definitions(convertAnimation_gtest PRIVATE
    PGO_TEST_CONVERT_ANIMATION_BIN="$<TARGET_FILE:convertAnimation>")
  set_property(TARGET convertAnimation_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(convertAnimation_gtest)
endif()

if(TARGET mergeCloseVertices)
  add_executable(mergeCloseVertices_gtest mergeCloseVertices_gtest.cpp)
  target_link_libraries(mergeCloseVertices_gtest PRIVATE GTest::gtest_main)
  add_dependencies(mergeCloseVertices_gtest mergeCloseVertices)
  target_compile_definitions(mergeCloseVertices_gtest PRIVATE
    PGO_TEST_MERGE_CLOSE_VERTICES_BIN="$<TARGET_FILE:mergeCloseVertices>")
  set_property(TARGET mergeCloseVertices_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(mergeCloseVertices_gtest)
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

if(PGO_ENABLE_CUDA)
  add_executable(cudaElasticEquivalence cudaElasticEquivalence.cpp)
  target_link_libraries(cudaElasticEquivalence PRIVATE simulationRunner cudaSolidDeformationModel cudaKoiterDeformationModel CUDA::cudart)
  target_include_directories(cudaElasticEquivalence PRIVATE ${CMAKE_SOURCE_DIR}/tests/src/core/solidDeformationModel)
  set_property(TARGET cudaElasticEquivalence PROPERTY FOLDER "tests/measurements")

  add_executable(cudaIPCSim_gtest cudaIPCSim_gtest.cpp)
  target_link_libraries(cudaIPCSim_gtest PRIVATE GTest::gtest_main simulationRunner CUDA::cudart)
  add_dependencies(cudaIPCSim_gtest validateIPCTrajectory libpgo_test_cubic_assets)
  target_compile_definitions(cudaIPCSim_gtest PRIVATE
    PGO_TEST_VALIDATE_IPC_TRAJECTORY_BIN="$<TARGET_FILE:validateIPCTrajectory>"
    PGO_CUDA_TEST_BOX_DIR="${CMAKE_SOURCE_DIR}/examples/assets/volume/box"
    PGO_CUDA_TEST_CUBIC_BOX_VEG="${LIBPGO_TEST_CUBIC_BOX_VEG}"
    PGO_CUDA_TEST_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/tests/assets/cudaIPC")
  set_property(TARGET cudaIPCSim_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(cudaIPCSim_gtest PROPERTIES LABELS "cuda")
endif()
