set(SIMULATION_RUNNER_SOURCES
  simulationRunner.cpp
  sampledSimulationRunner.cpp
  ipcSimulationRunner.cpp
  runSimCliLogging.cpp
  runSimFEMSetup.cpp
  runSimVolumeMeshIO.cpp
  runIPCExternalObjects.cpp
  runSimVertexTrajectory.cpp
  runIPCSimSetup.cpp)

set(SIMULATION_RUNNER_HEADERS
  simulationRunner.h
  runSimCliLogging.h
  runSimTiming.h
  runSimFEMSetup.h
  runSimVolumeMeshIO.h
  runIPCExternalObjects.h
  runSimVertexTrajectory.h
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

if(PGO_ENABLE_CUDA)
  target_sources(simulationRunner PRIVATE cudaIpcSimulationRunner.cpp cudaIpcSimulationRunner.h)
  target_compile_definitions(simulationRunner PRIVATE PGO_HAS_CUDA_IPC)
  target_link_libraries(simulationRunner PRIVATE cudaNonlinearOptimization cudaSolidDeformationModel cudaKoiterDeformationModel cudaContact)
endif()
