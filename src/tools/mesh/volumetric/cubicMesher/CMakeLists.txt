set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_library(cubicMesherCore STATIC
  cubicMesherIO.cpp
  triangleMeshVoxelizer.cpp)
target_include_directories(cubicMesherCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cubicMesherCore PUBLIC ${AVAILABLE_LIBS})
target_link_libraries(cubicMesherCore PRIVATE compilation_flag)

if(PGO_RELEASE_MODE_DEBUG)
  target_link_libraries(cubicMesherCore PRIVATE compilation_flag_for_debug)
endif()

set_property(TARGET cubicMesherCore PROPERTY FOLDER "tools/internal")

add_libpgo_tools(cubicMesher cubicMesher.cpp "")
target_link_libraries(cubicMesher PRIVATE cubicMesherCore argparse::argparse)
