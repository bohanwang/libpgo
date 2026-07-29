add_executable(configFileJSON_gtest configFileJSON_gtest.cpp)
target_link_libraries(configFileJSON_gtest PRIVATE GTest::gtest_main configFileJSON)
target_compile_definitions(configFileJSON_gtest PRIVATE
  LIBPGO_TEST_BOX_JSON="${LIBPGO_TEST_TET_BOX_JSON}"
  LIBPGO_TEST_CUBIC_BOX_JSON="${LIBPGO_TEST_CUBIC_BOX_JSON}"
  LIBPGO_TEST_SHELL_JSON="${CMAKE_SOURCE_DIR}/examples/configs/shell/shell-dynamic-sampled.json")
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

add_subdirectory(basicAlgorithms)
add_subdirectory(contact)
add_subdirectory(solidDeformationModel)
