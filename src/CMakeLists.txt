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
  animationIO
  CACHE INTERNAL ""
)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  add_subdirectory(core/${lib})
endforeach()

if(PGO_UI)
  add_subdirectory(render)
endif()

add_subdirectory(c)

if(PGO_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

add_subdirectory(tools)

add_subdirectory(tests)
