project(macros LANGUAGES CXX C)

add_library(macros INTERFACE)
target_link_libraries(macros INTERFACE compilation_flag)
target_include_directories(macros INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


