add_executable(configFileJSON_gtest configFileJSON_gtest.cpp)
target_link_libraries(configFileJSON_gtest PRIVATE GTest::gtest_main configFileJSON)
target_compile_definitions(configFileJSON_gtest PRIVATE
  LIBPGO_TEST_BOX_JSON="${LIBPGO_TEST_TET_BOX_JSON}"
  LIBPGO_TEST_CUBIC_BOX_JSON="${LIBPGO_TEST_CUBIC_BOX_JSON}"
  LIBPGO_TEST_SHELL_JSON="${CMAKE_SOURCE_DIR}/examples/configs/shell/shell_hang/shell-dynamic-ipc.json")
set_property(TARGET configFileJSON_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(configFileJSON_gtest)

add_executable(newtonSolver_gtest newtonSolver_gtest.cpp)
target_link_libraries(newtonSolver_gtest PRIVATE GTest::gtest_main nonlinearOptimization genericPotentialEnergies solidDeformationModel contact)
set_property(TARGET newtonSolver_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(newtonSolver_gtest)

if(APPLE AND TARGET BLAS::BLAS)
  add_executable(accelerateRuntimeSmoke accelerateRuntimeSmoke.c)
  target_link_libraries(accelerateRuntimeSmoke PRIVATE BLAS::BLAS)
  set_property(TARGET accelerateRuntimeSmoke PROPERTY FOLDER "tests/smoke")
  add_test(NAME accelerateRuntimeSmoke COMMAND accelerateRuntimeSmoke)
endif()

if (TARGET animationIO)
  add_executable(animationLoader_gtest animationLoader_gtest.cpp)
  target_link_libraries(animationLoader_gtest PRIVATE GTest::gtest_main animationIO)
  set_property(TARGET animationLoader_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(animationLoader_gtest)
endif()

add_executable(constraintFunctionsAssembler_gtest constraintFunctionsAssembler_gtest.cpp)
target_link_libraries(constraintFunctionsAssembler_gtest PRIVATE GTest::gtest_main nonlinearOptimization)
set_property(TARGET constraintFunctionsAssembler_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(constraintFunctionsAssembler_gtest)

add_executable(basicIO_gtest basicIO_gtest.cpp)
target_link_libraries(basicIO_gtest PRIVATE GTest::gtest_main basicIO)
set_property(TARGET basicIO_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(basicIO_gtest)

add_executable(eigenSupportNPY_gtest eigenSupportNPY_gtest.cpp)
target_link_libraries(eigenSupportNPY_gtest PRIVATE GTest::gtest_main eigenSupport)
set_property(TARGET eigenSupportNPY_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(eigenSupportNPY_gtest)

if(TARGET geometryPotentialEnergies)
  add_executable(centerOfMassMatchingEnergy_gtest centerOfMassMatchingEnergy_gtest.cpp)
  target_link_libraries(centerOfMassMatchingEnergy_gtest PRIVATE GTest::gtest_main geometryPotentialEnergies nonlinearOptimization mesh)
  set_property(TARGET centerOfMassMatchingEnergy_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(centerOfMassMatchingEnergy_gtest)
endif()

add_executable(barycentricCoordinates_gtest barycentricCoordinates_gtest.cpp)
target_link_libraries(barycentricCoordinates_gtest PRIVATE GTest::gtest_main interpolationCoordinates volumetricMesh)
set_property(TARGET barycentricCoordinates_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(barycentricCoordinates_gtest)

if(TARGET Gmsh::Gmsh)
  add_executable(loadMshFile_gtest loadMshFile_gtest.cpp)
  target_link_libraries(loadMshFile_gtest PRIVATE GTest::gtest_main volumetricMesh Gmsh::Gmsh)
  set_property(TARGET loadMshFile_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(loadMshFile_gtest)
endif()

if(TARGET tetgenInterface)
  add_executable(tetgenInterface_gtest tetgenInterface_gtest.cpp)
  target_link_libraries(tetgenInterface_gtest PRIVATE GTest::gtest_main tetgenInterface mesh)
  set_property(TARGET tetgenInterface_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(tetgenInterface_gtest)
endif()

add_subdirectory(basicAlgorithms)
add_subdirectory(contact)
add_subdirectory(solidDeformationModel)
