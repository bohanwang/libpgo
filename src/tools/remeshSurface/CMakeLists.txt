set(AVAILABLE_LIBS)

if(NOT TARGET cgalInterface OR NOT TARGET geogramInterface)
  return()
endif()

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_libpgo_tools(remeshSurface remeshSurface.cpp "${AVAILABLE_LIBS}")

target_link_libraries(remeshSurface PRIVATE argparse::argparse)