set(SIMULATION_RUNNER_SOURCES
  simulationRunner.cpp
  sampledSimulationRunner.cpp
  ipcSimulationRunner.cpp
  runSimCliLogging.cpp
  runSimFEMSetup.cpp
  runSimVolumeMeshIO.cpp
  runIPCSimSetup.cpp)

set(SIMULATION_RUNNER_HEADERS
  simulationRunner.h
  runSimCliLogging.h
  runSimFEMSetup.h
  runSimVolumeMeshIO.h
  runIPCSimSetup.h)

add_libpgo_lib(simulationRunner "${SIMULATION_RUNNER_SOURCES}" "${SIMULATION_RUNNER_HEADERS}")

target_include_directories(simulationRunner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(simulationRunner PUBLIC
  basicIO
  configFileJSON
  constraintPotentialEnergies
  contact
  eigenSupport
  genericPotentialEnergies
  interpolationCoordinates
  libiglInterface
  mesh
  nonlinearOptimization
  pgoLogging
  simulation
  solidDeformationModel
  volumetricMesh)
