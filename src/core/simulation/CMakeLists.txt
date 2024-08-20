set(SIMULATION_HEADERS
  timeIntegrator.h
  timeIntegratorSolver.h
  implicitBackwardEulerTimeIntegratorHelper.h
  implicitBackwardEulerTimeIntegrator.h
  TRBDF2TimeIntegratorHelper.h
  TRBDF2TimeIntegrator.h
)

set(SIMULATION_SOURCES
  timeIntegrator.cpp
  timeIntegratorSolver.cpp
  implicitBackwardEulerTimeIntegratorHelper.cpp
  implicitBackwardEulerTimeIntegrator.cpp
  TRBDF2TimeIntegratorHelper.cpp
  TRBDF2TimeIntegrator.cpp
)

set(SIMULATION_DEPS constraintPotentialEnergies solidDeformationModel)

foreach(tgt ${SIMULATION_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "simulation is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(simulation "${SIMULATION_SOURCES}" "${SIMULATION_HEADERS}")

target_link_libraries(simulation PUBLIC ${SIMULATION_DEPS})