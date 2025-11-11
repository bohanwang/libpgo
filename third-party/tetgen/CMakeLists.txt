# Add a library to the project using the specified source files.
# In Linux/Unix, it will creates the libtet.a
add_library(tetgen STATIC tetgen.h tetgen.cxx predicates.cxx)

if (TARGET CGAL::CGAL)
  target_compile_definitions(tetgen PRIVATE USE_CGAL_PREDICATES)
  target_link_libraries(tetgen PUBLIC CGAL::CGAL)
endif()

target_compile_definitions(tetgen PUBLIC TETLIBRARY)

target_include_directories(tetgen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


set_property(TARGET tetgen PROPERTY FOLDER "third-party")