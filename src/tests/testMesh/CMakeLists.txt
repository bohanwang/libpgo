set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()


add_libpgo_tools(testMesh testMesh.cpp "${AVAILABLE_LIBS}")

set_property(TARGET testMesh PROPERTY FOLDER "tests")