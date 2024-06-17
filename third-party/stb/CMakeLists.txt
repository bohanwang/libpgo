project(stb LANGUAGES CXX)

find_path(STB_INCLUDE_DIRS "stb_image.h")

add_library(stb STATIC stb_libs.cpp)
target_include_directories(stb PRIVATE ${STB_INCLUDE_DIRS})


set_property(TARGET stb PROPERTY FOLDER "3rdParty")