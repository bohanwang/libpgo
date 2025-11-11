set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_libpgo_tools(tetMesher tetMesher.cpp "${AVAILABLE_LIBS}")
target_link_libraries(tetMesher PRIVATE argparse::argparse)

add_libpgo_tools(generateTetMeshSurfaceMesh generateTetMeshSurfaceMesh.cpp "${AVAILABLE_LIBS}")

if(TARGET Gmsh::Gmsh)
  add_libpgo_tools(mshFileToVegFile mshFileToVegFile.cpp "${AVAILABLE_LIBS}")
endif()
