cmake_minimum_required(VERSION 3.25.0 FATAL_ERROR)

project(pgo_c_test C CXX)

find_package(pgo CONFIG REQUIRED)


add_executable(pgo_c_test pgo_c_test.cpp)

target_link_libraries(pgo_c_test PRIVATE pgo::pgo_c)