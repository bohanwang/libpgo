set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

if(NOT TARGET libiglInterface)
  return()
endif()

add_libpgo_tools(generateFBMSUnionSurface generateFBMSUnionSurface.cpp "${AVAILABLE_LIBS}")
target_link_libraries(generateFBMSUnionSurface PRIVATE argparse::argparse)
