set(CONTACT_HEADERS
  CCDKernel.h
  contactEnergyUtilities.h
  pointPenetrationBarrierEnergy.h
  pointPenetrationEnergy.h
  pointTrianglePairBarrierEnergy.h
  pointTrianglePairCouplingEnergyWithCollision.h
  triangleMeshExternalContactHandler.h
  triangleMeshSelfContactDetection.h
  triangleMeshSelfContactHandler.h
)

set(CONTACT_SOURCES
  CCDKernel.cpp
  contactEnergyUtilities.cpp
  pointPenetrationBarrierEnergy.cpp
  pointPenetrationEnergy.cpp
  pointTrianglePairBarrierEnergy.cpp
  pointTrianglePairCouplingEnergyWithCollision.cpp
  triangleMeshExternalContactHandler.cpp
  triangleMeshSelfContactDetection.cpp
  triangleMeshSelfContactHandler.cpp
)

set(CONTACT_DEPS mesh nonlinearOptimization CCD::SafeCCD CCD::Exact)

foreach(tgt ${CONTACT_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "selfContact is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(contact "${CONTACT_SOURCES}" "${CONTACT_HEADERS}")

target_link_libraries(contact PUBLIC ${CONTACT_DEPS})
