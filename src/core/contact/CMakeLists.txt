set(CONTACT_HEADERS
  CCDKernel.h
  contactEnergyUtilities.h
  mappedSurfacePotentialEnergy.h
  embeddedSurfaceFloorPotentialEnergy.h
  embeddedSurfaceIPCPotentialEnergy.h
  ipc/broadPhase/spatialHashGrid.h
  ipc/broadPhase/surfaceIPCSelfBroadPhase.h
  ipc/core/surfaceIPCBarrierAssembler.h
  ipc/core/surfaceIPCMaxStep.h
  ipc/core/surfaceIPCCore.h
  ipc/core/surfaceIPCPairs.h
  ipc/geometry/ipcBarrier.h
  ipc/geometry/ipcCCD.h
  ipc/geometry/ipcDistancePrimitives.h
  ipc/geometry/ipcHessianProjection.h
  ipc/geometry/generated/CIPC_autogen.h
  ipc/geometry/generated/CIPC_autogen_ll.h
  ipc/profiling/surfaceIPCProfiling.h
  ipc/topology/surfaceIPCTopology.h
  pointPenetrationEnergy.h
  pointTrianglePairCouplingEnergyWithCollision.h
  triangleMeshExternalContactHandler.h
  triangleMeshSelfContactDetection.h
  triangleMeshSelfContactHandler.h
)

set(CONTACT_SOURCES
  CCDKernel.cpp
  contactEnergyUtilities.cpp
  mappedSurfacePotentialEnergy.cpp
  embeddedSurfaceFloorPotentialEnergy.cpp
  embeddedSurfaceIPCPotentialEnergy.cpp
  ipc/broadPhase/spatialHashGrid.cpp
  ipc/broadPhase/surfaceIPCSelfBroadPhase.cpp
  ipc/core/surfaceIPCBarrierAssembler.cpp
  ipc/core/surfaceIPCMaxStep.cpp
  ipc/core/surfaceIPCCore.cpp
  ipc/geometry/ipcBarrier.cpp
  ipc/geometry/ipcCCD.cpp
  ipc/geometry/ipcDistancePrimitives.cpp
  ipc/geometry/ipcHessianProjection.cpp
  ipc/topology/surfaceIPCTopology.cpp
  pointPenetrationEnergy.cpp
  pointTrianglePairCouplingEnergyWithCollision.cpp
  triangleMeshExternalContactHandler.cpp
  triangleMeshSelfContactDetection.cpp
  triangleMeshSelfContactHandler.cpp
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
