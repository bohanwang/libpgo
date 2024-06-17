project(glad LANGUAGES C)

add_library(glad STATIC ./src/glad.c)
target_link_libraries(glad PRIVATE compilation_flag)
target_include_directories(glad PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

set_property(TARGET glad PROPERTY FOLDER "3rdParty")