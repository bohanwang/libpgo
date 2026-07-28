set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_libpgo_tools(runSim "runSim.cpp" "simulationRunner")

target_link_libraries(runSim PRIVATE argparse::argparse)

add_libpgo_tools(runShellSim "runShellSim.cpp;${CMAKE_SOURCE_DIR}/src/simulationRunner/runSimCliLogging.cpp" "${AVAILABLE_LIBS}")

target_link_libraries(runShellSim PRIVATE argparse::argparse)
target_include_directories(runShellSim PRIVATE ${CMAKE_SOURCE_DIR}/src/simulationRunner)

add_libpgo_tools(runIPCSim "runIPCSim.cpp" "simulationRunner")

target_link_libraries(runIPCSim PRIVATE argparse::argparse)
