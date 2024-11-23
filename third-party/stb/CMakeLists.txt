add_library(stb STATIC stb_libs.cpp stb_image.h stb_image_write.h)
target_include_directories(stb PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set_property(TARGET stb PROPERTY FOLDER "third-party")
