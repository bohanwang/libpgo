project(tiny_obj_loader LANGUAGES C)

add_library(tiny_obj_loader STATIC tiny_obj_loader.cc tiny_obj_loader.h)
target_link_libraries(tiny_obj_loader PRIVATE compilation_flag)

target_include_directories(tiny_obj_loader PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(tiny_obj_loader PUBLIC TINYOBJLOADER_USE_DOUBLE)

set_property(TARGET tiny_obj_loader PROPERTY FOLDER "third-party")
