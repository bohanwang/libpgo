add_executable(configFileJSON_gtest configFileJSON_gtest.cpp)
target_link_libraries(configFileJSON_gtest PRIVATE GTest::gtest_main configFileJSON)
target_compile_definitions(configFileJSON_gtest PRIVATE
  LIBPGO_TEST_BOX_JSON="${CMAKE_SOURCE_DIR}/examples/box/box.json"
  LIBPGO_TEST_CUBIC_BOX_JSON="${CMAKE_SOURCE_DIR}/examples/cubic/box/box.json"
  LIBPGO_TEST_SHELL_JSON="${CMAKE_SOURCE_DIR}/examples/shell/shell.json")
set_property(TARGET configFileJSON_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(configFileJSON_gtest)

if (TARGET animationIO)
  add_executable(animationLoader_gtest animationLoader_gtest.cpp)
  target_link_libraries(animationLoader_gtest PRIVATE GTest::gtest_main animationIO)
  set_property(TARGET animationLoader_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(animationLoader_gtest)
endif()

if(TARGET nonlinearOptimization)
  add_executable(NewtonSolver_gtest NewtonSolver_gtest.cpp)
  target_link_libraries(NewtonSolver_gtest PRIVATE GTest::gtest_main nonlinearOptimization)
  set_property(TARGET NewtonSolver_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(NewtonSolver_gtest)
endif()

if(TARGET simulation)
  add_executable(implicitBackwardEulerTimeIntegrator_gtest implicitBackwardEulerTimeIntegrator_gtest.cpp)
  target_link_libraries(implicitBackwardEulerTimeIntegrator_gtest PRIVATE GTest::gtest_main simulation)
  set_property(TARGET implicitBackwardEulerTimeIntegrator_gtest PROPERTY FOLDER "tests/gtest")
  gtest_discover_tests(implicitBackwardEulerTimeIntegrator_gtest)
endif()

add_subdirectory(basicAlgorithms)
add_subdirectory(contact)
add_subdirectory(profiling)
add_subdirectory(solidDeformationModel)
