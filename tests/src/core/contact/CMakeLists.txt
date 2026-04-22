add_executable(cipcProfiling_gtest cipcProfiling_gtest.cpp)
target_link_libraries(cipcProfiling_gtest PRIVATE GTest::gtest_main contact profiling)
set_property(TARGET cipcProfiling_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(cipcProfiling_gtest)
