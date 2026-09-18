add_libpgo_tools(runSim "runSim.cpp" "simulationRunner")

target_link_libraries(runSim PRIVATE argparse::argparse)

add_libpgo_tools(runIPCSim "runIPCSim.cpp" "simulationRunner")

target_link_libraries(runIPCSim PRIVATE argparse::argparse)
