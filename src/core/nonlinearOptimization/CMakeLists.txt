set(NONLINEAROPTIMIZATION_HEADERS 
  automaticDifferentiation_autodiff.h

  potentialEnergy.h 
  potentialEnergies.h 

  potentialEnergyFromConstraintFunctions.h
  
  constraintFunction.h constraintFunctions.h constraintFunctionsAssember.h 
  linearConstraintFunctions.h
  
  nonlinearProblem.h lagrangian.h 
  NewtonRaphsonSolver.h minimizeEnergy.h
  finiteDifference.h
  determinantDerivatives.h polarDecompositionDerivatives.h naturalCubicSplineDerivatives.h svdDerivatives.h 
  barrierFunction.h
  naturalCubicSplineFitting.h
  lineSearch.h
)

set(NONLINEAROPTIMIZATION_SOURCES 
  potentialEnergies.cpp
  
  constraintFunctions.cpp constraintFunctionsAssember.cpp 
  linearConstraintFunctions.cpp

  potentialEnergyFromConstraintFunctions.cpp
  nonlinearProblem.cpp lagrangian.cpp  
  NewtonRaphsonSolver.cpp minimizeEnergy.cpp
  finiteDifference.cpp    
  determinantDerivatives.cpp polarDecompositionDerivatives.cpp naturalCubicSplineDerivatives.cpp svdDerivatives.cpp
  barrierFunction.cpp
  naturalCubicSplineFitting.cpp
  lineSearch.cpp
)

if(TARGET Knitro::Knitro)
  set(NONLINEAROPTIMIZATION_HEADERS ${NONLINEAROPTIMIZATION_HEADERS} knitroProblem.h knitroOptimizer.h)
  set(NONLINEAROPTIMIZATION_SOURCES ${NONLINEAROPTIMIZATION_SOURCES} knitroProblem.cpp knitroOptimizer.cpp)
endif()

if(TARGET Ipopt::Ipopt)
  set(NONLINEAROPTIMIZATION_HEADERS ${NONLINEAROPTIMIZATION_HEADERS} IpoptProblem.h IpoptOptimizer.h)
  set(NONLINEAROPTIMIZATION_SOURCES ${NONLINEAROPTIMIZATION_SOURCES} IpoptProblem.cpp IpoptOptimizer.cpp)
endif()

set(NONLINEAROPTIMIZATION_DEPS eigenSupport pgoLogging autodiff::autodiff)

foreach(tgt ${NONLINEAROPTIMIZATION_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "nonlinearOptimization is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(nonlinearOptimization "${NONLINEAROPTIMIZATION_SOURCES}" "${NONLINEAROPTIMIZATION_HEADERS}")

foreach(lib ${NONLINEAROPTIMIZATION_DEPS})
  target_link_libraries(nonlinearOptimization PUBLIC ${lib})
endforeach()

if(TARGET Knitro::Knitro)
  target_link_libraries(nonlinearOptimization PRIVATE Knitro::Knitro)
endif()

if(TARGET Ipopt::Ipopt)
  target_link_libraries(nonlinearOptimization PRIVATE Ipopt::Ipopt)
endif()
