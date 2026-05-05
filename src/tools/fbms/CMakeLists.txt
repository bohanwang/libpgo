set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

if(TARGET pgoLogging)
  add_libpgo_tools(computeStressFieldStats computeStressFieldStats.cpp "pgoLogging")
  target_link_libraries(computeStressFieldStats PRIVATE argparse::argparse nlohmann_json::nlohmann_json fmt::fmt)
endif()

if(NOT TARGET libiglInterface)
  return()
endif()

add_libpgo_tools(generateFBMSUnionSurface generateFBMSUnionSurface.cpp "${AVAILABLE_LIBS}")
target_link_libraries(generateFBMSUnionSurface PRIVATE argparse::argparse)
