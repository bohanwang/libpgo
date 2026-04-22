set(CONTACT_HEADERS
  CCDKernel.h
  contactEnergyUtilities.h
  surfaceIPCProfiling.h
  pointPenetrationEnergy.h
  pointTrianglePairCouplingEnergyWithCollision.h
  triangleMeshExternalContactHandler.h
  triangleMeshSelfContactDetection.h
  triangleMeshSelfContactHandler.h
  CIPC.h
  CIPC_autogen.h
  CIPC_autogen_ll.h
)

set(CONTACT_SOURCES
  CCDKernel.cpp
  contactEnergyUtilities.cpp
  pointPenetrationEnergy.cpp
  pointTrianglePairCouplingEnergyWithCollision.cpp
  triangleMeshExternalContactHandler.cpp
  triangleMeshSelfContactDetection.cpp
  triangleMeshSelfContactHandler.cpp
  CIPC.cpp
)

set(CONTACT_DEPS mesh nonlinearOptimization profiling CCD::SafeCCD CCD::Exact)

foreach(tgt ${CONTACT_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "selfContact is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(contact "${CONTACT_SOURCES}" "${CONTACT_HEADERS}")

target_link_libraries(contact PUBLIC ${CONTACT_DEPS})
