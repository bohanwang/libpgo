set(SIMULATION_API_HEADERS
  runSimCore.h
)

set(SIMULATION_API_SOURCES
  runSimCore.cpp
)

set(SIMULATION_API_DEPS
  basicIO
  configFileJSON
  contact
  constraintPotentialEnergies
  fileService
  genericPotentialEnergies
  interpolationCoordinates
  mesh
  nonlinearOptimization
  pgoLogging
  simulation
  solidDeformationModel
  volumetricMesh
)

foreach(tgt ${SIMULATION_API_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "simulationApi is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(simulationApi "${SIMULATION_API_SOURCES}" "${SIMULATION_API_HEADERS}")
target_link_libraries(simulationApi PUBLIC ${SIMULATION_API_DEPS})

add_subdirectory(c)
if(PGO_FEATURE_PYTHON)
  add_subdirectory(python)
endif()
