add_executable(cubicMesher_gtest cubicMesher_gtest.cpp)
target_link_libraries(cubicMesher_gtest PRIVATE GTest::gtest_main cubicMesherCore)
target_include_directories(cubicMesher_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/tools/cubicMesher)
add_dependencies(cubicMesher_gtest cubicMesher)
target_compile_definitions(cubicMesher_gtest PRIVATE
  PGO_TEST_CUBIC_MESHER_BIN="$<TARGET_FILE:cubicMesher>")
set_property(TARGET cubicMesher_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cubicMesher_gtest)

add_executable(tetMesher_gtest tetMesher_gtest.cpp)
target_link_libraries(tetMesher_gtest PRIVATE GTest::gtest_main volumetricMesh mesh)
add_dependencies(tetMesher_gtest tetMesher)
target_compile_definitions(tetMesher_gtest PRIVATE
  PGO_TEST_TET_MESHER_BIN="$<TARGET_FILE:tetMesher>"
  PGO_TEST_TET_MESHER_HAS_TET_WILD=$<BOOL:${PGO_TET_MESHER_USE_TET_WILD}>)
set_property(TARGET tetMesher_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(tetMesher_gtest)

add_executable(meshQualityCheck_gtest meshQualityCheck_gtest.cpp)
target_link_libraries(meshQualityCheck_gtest PRIVATE GTest::gtest_main nlohmann_json::nlohmann_json)
add_dependencies(meshQualityCheck_gtest meshQualityCheck)
target_compile_definitions(meshQualityCheck_gtest PRIVATE
  PGO_TEST_MESH_QUALITY_CHECK_BIN="$<TARGET_FILE:meshQualityCheck>")
set_property(TARGET meshQualityCheck_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(meshQualityCheck_gtest)

if(TARGET convertAnimation)
  add_executable(convertAnimation_gtest convertAnimation_gtest.cpp)
  target_link_libraries(convertAnimation_gtest PRIVATE GTest::gtest_main)
  add_dependencies(convertAnimation_gtest convertAnimation)
  target_compile_definitions(convertAnimation_gtest PRIVATE
    PGO_TEST_CONVERT_ANIMATION_BIN="$<TARGET_FILE:convertAnimation>")
  set_property(TARGET convertAnimation_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(convertAnimation_gtest)
endif()

add_executable(runSimShared_gtest
  runSimShared_gtest.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimCliLogging.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimVolumeMeshIO.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimFEMSetup.cpp)
target_link_libraries(runSimShared_gtest PRIVATE
  GTest::gtest_main
  configFileJSON
  contact
  interpolationCoordinates
  solidDeformationModel
  volumetricMesh
  mesh)
target_include_directories(runSimShared_gtest PRIVATE
  ${CMAKE_SOURCE_DIR}/src/tools/runSim
  ${CMAKE_SOURCE_DIR}/src/core/solidDeformationModel)
target_compile_definitions(runSimShared_gtest PRIVATE
  LIBPGO_TEST_TET_BOX_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/tet/box/box.veg"
  LIBPGO_TEST_TET_BOX_OBJ="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/tet/box/box.obj"
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/cubic/box/box.veg"
  LIBPGO_TEST_CUBIC_BOX_OBJ="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/cubic/box/box.obj")
set_property(TARGET runSimShared_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runSimShared_gtest)

add_executable(runIPCSim_gtest
  runIPCSim_gtest.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runIPCSimCli.cpp)
target_link_libraries(runIPCSim_gtest PRIVATE
  GTest::gtest_main
  runIPCSimCore
  argparse::argparse
  configFileJSON
  contact
  constraintPotentialEnergies
  interpolationCoordinates
  libiglInterface
  simulation
  solidDeformationModel
  volumetricMesh
  mesh)
target_include_directories(runIPCSim_gtest PRIVATE
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
  LIBPGO_TEST_SHELL_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures/shell"
  LIBPGO_TEST_LEGACY_TET_BOX_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/tet/box"
  LIBPGO_TEST_LEGACY_CUBIC_BOX_DIR="${CMAKE_SOURCE_DIR}/tests/fixtures/legacy/cubic/box"
  LIBPGO_TEST_IPC_TET_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/examples/ipc/tet/box-hang"
  LIBPGO_TEST_IPC_CUBIC_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/examples/ipc/cubic/box-hang"
  PGO_TEST_RUN_IPC_SIM_BIN="$<TARGET_FILE:runIPCSim>")
add_dependencies(runIPCSim_gtest runIPCSim)
set_property(TARGET runIPCSim_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runIPCSim_gtest)

if(TARGET rawSurfaceCleanup)
  add_executable(rawSurfaceCleanup_gtest rawSurfaceCleanup_gtest.cpp)
  target_link_libraries(rawSurfaceCleanup_gtest PRIVATE GTest::gtest_main nlohmann_json::nlohmann_json)
  add_dependencies(rawSurfaceCleanup_gtest rawSurfaceCleanup)
  target_compile_definitions(rawSurfaceCleanup_gtest PRIVATE
    PGO_TEST_RAW_SURFACE_CLEANUP_BIN="$<TARGET_FILE:rawSurfaceCleanup>")
  set_property(TARGET rawSurfaceCleanup_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(rawSurfaceCleanup_gtest)
endif()

if(TARGET remeshSurface AND TARGET meshQualityCheck)
  add_executable(remeshSurface_gtest remeshSurface_gtest.cpp)
  target_link_libraries(remeshSurface_gtest PRIVATE GTest::gtest_main nlohmann_json::nlohmann_json)
  add_dependencies(remeshSurface_gtest remeshSurface meshQualityCheck)
  target_compile_definitions(remeshSurface_gtest PRIVATE
    PGO_TEST_REMESH_SURFACE_BIN="$<TARGET_FILE:remeshSurface>"
    PGO_TEST_MESH_QUALITY_CHECK_BIN="$<TARGET_FILE:meshQualityCheck>")
  set_property(TARGET remeshSurface_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(remeshSurface_gtest)
endif()
