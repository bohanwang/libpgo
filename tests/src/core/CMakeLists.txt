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

add_subdirectory(solidDeformationModel)
