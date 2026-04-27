set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

set(TET_MESHER_SOURCES
  tetMesher.cpp
  tetMesherBackend.cpp
  tetgenBackend.cpp)

if(PGO_TET_MESHER_USE_TET_WILD)
  list(APPEND TET_MESHER_SOURCES tetwildBackend.cpp)
else()
  list(APPEND TET_MESHER_SOURCES tetwildBackendDisabled.cpp)
endif()

add_libpgo_tools(tetMesher "${TET_MESHER_SOURCES}" "${AVAILABLE_LIBS}")
target_link_libraries(tetMesher PRIVATE argparse::argparse configFileJSON)

if(PGO_TET_MESHER_USE_TET_WILD)
  target_link_libraries(tetMesher PRIVATE FloatTetwild)
  target_compile_definitions(tetMesher PRIVATE PGO_TET_MESHER_HAS_TET_WILD=1)
endif()

add_libpgo_tools(generateTetMeshSurfaceMesh generateTetMeshSurfaceMesh.cpp "${AVAILABLE_LIBS}")

if(TARGET Gmsh::Gmsh)
  add_libpgo_tools(mshFileToVegFile mshFileToVegFile.cpp "${AVAILABLE_LIBS}")
endif()
