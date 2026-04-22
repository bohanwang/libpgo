set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_libpgo_tools(runSim "runSim.cpp;runSimVolumeMeshIO.cpp;runSimFEMSetup.cpp;runSimCliLogging.cpp" "${AVAILABLE_LIBS}")

target_link_libraries(runSim PRIVATE argparse::argparse)

add_libpgo_tools(runShellSim "runShellSim.cpp;runSimCliLogging.cpp" "${AVAILABLE_LIBS}")

target_link_libraries(runShellSim PRIVATE argparse::argparse)

add_libpgo_tools(runIPCSim "runIPCSim.cpp;runSimCliLogging.cpp" "${AVAILABLE_LIBS}")

target_link_libraries(runIPCSim PRIVATE argparse::argparse)
