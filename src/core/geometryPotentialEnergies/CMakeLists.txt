set(GEOMETRY_POTENTIAL_ENERGIES_HEADERS
  surfaceSmoothnessAbsoluteMeanCurvature.h
  surfaceTriangleDeformation.h
  vertexAffineToPositionEnergy.h
  triangleAffineToPositionEnergy.h
  centerOfMassMatchingEnergy.h
)

set(GEOMETRY_POTENTIAL_ENERGIES_SOURCES
  surfaceSmoothnessAbsoluteMeanCurvature.cpp
  surfaceTriangleDeformation.cpp
  vertexAffineToPositionEnergy.cpp
  triangleAffineToPositionEnergy.cpp
  centerOfMassMatchingEnergy.cpp
)

set(GEOMETRY_POTENTIAL_ENERGIES_DEPS nonlinearOptimization mesh)

if(TARGET MKL::MKL)
  list(APPEND GEOMETRY_POTENTIAL_ENERGIES_HEADERS smoothRSEnergy.h)
  list(APPEND GEOMETRY_POTENTIAL_ENERGIES_SOURCES smoothRSEnergy.cpp)
  list(APPEND GEOMETRY_POTENTIAL_ENERGIES_DEPS MKL::MKL)
endif()

foreach(tgt ${GEOMETRY_POTENTIAL_ENERGIES_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "geometryPotentialEnergies is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(geometryPotentialEnergies "${GEOMETRY_POTENTIAL_ENERGIES_SOURCES}" "${GEOMETRY_POTENTIAL_ENERGIES_HEADERS}")

target_link_libraries(geometryPotentialEnergies PUBLIC ${GEOMETRY_POTENTIAL_ENERGIES_DEPS})