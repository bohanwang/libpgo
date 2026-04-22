add_executable(cipcProfiling_gtest cipcProfiling_gtest.cpp)
target_link_libraries(cipcProfiling_gtest PRIVATE GTest::gtest_main contact profiling)
set_property(TARGET cipcProfiling_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cipcProfiling_gtest)

add_executable(surfaceIPCCore_gtest surfaceIPCCore_gtest.cpp)
target_link_libraries(surfaceIPCCore_gtest PRIVATE GTest::gtest_main contact)
set_property(TARGET surfaceIPCCore_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(surfaceIPCCore_gtest)

add_executable(cipcPotentialEnergy_gtest cipcPotentialEnergy_gtest.cpp)
target_link_libraries(cipcPotentialEnergy_gtest PRIVATE GTest::gtest_main contact)
set_property(TARGET cipcPotentialEnergy_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cipcPotentialEnergy_gtest)
