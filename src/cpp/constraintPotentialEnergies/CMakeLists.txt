set(CONSTRAINT_POTENTIAL_ENERGIES_HEADERS
  potentialEnergyAligningMeshConnectivity.h
  multiVertexConstrainedRigidMotion.h
  multiVertexPullingSoftConstraints.h
  multiVertexSlidingSoftConstraints.h
  multiVertexPullingSoftConstraintsPOrder.h
  barycentricCoordinateSlidingSoftConstraints.h
)

set(CONSTRAINT_POTENTIAL_ENERGIES_SOURCES
  potentialEnergyAligningMeshConnectivity.cpp
  multiVertexConstrainedRigidMotion.cpp
  multiVertexPullingSoftConstraints.cpp
  multiVertexSlidingSoftConstraints.cpp
  multiVertexPullingSoftConstraintsPOrder.cpp
  barycentricCoordinateSlidingSoftConstraints.cpp
)

set(CONSTRAINT_POTENTIAL_ENERGIES_DEPS nonlinearOptimization mesh)

foreach(tgt ${CONSTRAINT_POTENTIAL_ENERGIES_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "constraintPotentialEnergies is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(constraintPotentialEnergies "${CONSTRAINT_POTENTIAL_ENERGIES_SOURCES}" "${CONSTRAINT_POTENTIAL_ENERGIES_HEADERS}")

target_link_libraries(constraintPotentialEnergies PUBLIC ${CONSTRAINT_POTENTIAL_ENERGIES_DEPS})