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
  geogramInterface
  libiglInterface
  tetgenInterface
  CACHE INTERNAL "Global library targets for the PGO project"
)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  add_subdirectory(core/${lib})
endforeach()

if(PGO_UI)
  message(FATAL_ERROR "Not supported yet")
endif()

add_subdirectory(c)

add_subdirectory(tools)

add_subdirectory(tests)
