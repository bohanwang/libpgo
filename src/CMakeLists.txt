set(PGO_GLOBAL_LIBRARY_TARGETS
  macros
  basicAlgorithms
  pgoLogging
  basicIO
  eigenSupport
  configFileJSON  
  mesh
  nonlinearOptimization
  genericPotentialEnergies
  geometryPotentialEnergies
  constraintPotentialEnergies
  contact
  volumetricMesh
  interpolationCoordinates
  solidDeformationModel
  simulation
  cgalInterface
)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  add_subdirectory(cpp/${lib})
endforeach()

add_subdirectory(c)

if(PGO_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

add_subdirectory(tools)

add_subdirectory(tests)
