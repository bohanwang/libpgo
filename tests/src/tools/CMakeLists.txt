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
  runSim_gtest.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimCliLogging.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimVolumeMeshIO.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimFEMSetup.cpp)
target_link_libraries(runSim_gtest PRIVATE
  GTest::gtest_main
  configFileJSON
  contact
  interpolationCoordinates
  solidDeformationModel
  volumetricMesh
  mesh)
target_include_directories(runSim_gtest PRIVATE
  ${CMAKE_SOURCE_DIR}/src/tools/runSim
  ${CMAKE_SOURCE_DIR}/src/core/solidDeformationModel)
target_compile_definitions(runSim_gtest PRIVATE
  LIBPGO_TEST_TET_BOX_VEG="${CMAKE_SOURCE_DIR}/examples/legacy/tet/box/box.veg"
  LIBPGO_TEST_TET_BOX_OBJ="${CMAKE_SOURCE_DIR}/examples/legacy/tet/box/box.obj"
  LIBPGO_TEST_CUBIC_BOX_VEG="${CMAKE_SOURCE_DIR}/examples/legacy/cubic/box/box.veg"
  LIBPGO_TEST_CUBIC_BOX_OBJ="${CMAKE_SOURCE_DIR}/examples/legacy/cubic/box/box.obj"
  LIBPGO_TEST_SHELL_JSON="${CMAKE_SOURCE_DIR}/examples/legacy/shell/shell.json"
  PGO_TEST_RUN_SHELL_SIM_BIN="$<TARGET_FILE:runShellSim>")
add_dependencies(runSim_gtest runShellSim)
set_property(TARGET runSim_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runSim_gtest)

add_executable(runIPCSim_gtest
  runIPCSim_gtest.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimCliLogging.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimVolumeMeshIO.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runSimFEMSetup.cpp
  ${CMAKE_SOURCE_DIR}/src/tools/runSim/runIPCSimSetup.cpp)
target_link_libraries(runIPCSim_gtest PRIVATE
  GTest::gtest_main
  configFileJSON
  contact
  constraintPotentialEnergies
  interpolationCoordinates
  libiglInterface
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
  LIBPGO_TEST_SHELL_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/examples/legacy/shell"
  LIBPGO_TEST_IPC_TET_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/examples/ipc/tet/box-hang"
  LIBPGO_TEST_IPC_CUBIC_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/examples/ipc/cubic/box-hang"
  PGO_TEST_RUN_IPC_SIM_BIN="$<TARGET_FILE:runIPCSim>")
add_dependencies(runIPCSim_gtest runIPCSim)
set_property(TARGET runIPCSim_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(runIPCSim_gtest)

if(TARGET generateFBMSUnionSurface)
  add_executable(generateFBMSUnionSurface_gtest generateFBMSUnionSurface_gtest.cpp)
  target_link_libraries(generateFBMSUnionSurface_gtest PRIVATE GTest::gtest_main mesh)
  add_dependencies(generateFBMSUnionSurface_gtest generateFBMSUnionSurface)
  target_compile_definitions(generateFBMSUnionSurface_gtest PRIVATE
    PGO_TEST_GENERATE_FBMS_UNION_SURFACE_BIN="$<TARGET_FILE:generateFBMSUnionSurface>")
  set_property(TARGET generateFBMSUnionSurface_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(generateFBMSUnionSurface_gtest)
endif()
