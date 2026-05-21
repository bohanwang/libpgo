set(CONTACT_HEADERS
  legacy_penalty/CCDKernel.h
  legacy_penalty/contactEnergyUtilities.h
  legacy_penalty/pointPenetrationEnergy.h
  legacy_penalty/pointTrianglePairCouplingEnergyWithCollision.h
  legacy_penalty/triangleMeshExternalContactHandler.h
  legacy_penalty/triangleMeshSelfContactDetection.h
  legacy_penalty/triangleMeshSelfContactHandler.h
  ipc/CIPC.h
  ipc/mappedSurfacePotentialEnergy.h
  ipc/embeddedSurfaceFloorPotentialEnergy.h
  ipc/embeddedSurfaceIPCPotentialEnergy.h
  ipc/broadPhase/spatialHashGrid.h
  ipc/broadPhase/surfaceIPCBroadPhase.h
  ipc/core/surfaceIPCSelfBarrierAssembler.h
  ipc/core/surfaceIPCExternalBarrierAssembler.h
  ipc/core/surfaceIPCBarrierKernels.h
  ipc/core/surfaceIPCMaxStep.h
  ipc/core/surfaceIPCCore.h
  ipc/core/surfaceIPCActiveSet.h
  ipc/core/surfaceIPCPairs.h
  ipc/external/obstacleSurface.h
  ipc/external/obstaclePoseCache.h
  ipc/geometry/ipcBarrier.h
  ipc/geometry/ipcCCD.h
  ipc/geometry/ipcDistancePrimitives.h
  ipc/geometry/ipcHessianProjection.h
  ipc/profiling/surfaceIPCProfiling.h
  ipc/topology/surfaceIPCTopology.h
)

set(CONTACT_SOURCES
  legacy_penalty/CCDKernel.cpp
  legacy_penalty/contactEnergyUtilities.cpp
  legacy_penalty/pointPenetrationEnergy.cpp
  legacy_penalty/pointTrianglePairCouplingEnergyWithCollision.cpp
  legacy_penalty/triangleMeshExternalContactHandler.cpp
  legacy_penalty/triangleMeshSelfContactDetection.cpp
  legacy_penalty/triangleMeshSelfContactHandler.cpp
  ipc/CIPC.cpp
  ipc/mappedSurfacePotentialEnergy.cpp
  ipc/embeddedSurfaceFloorPotentialEnergy.cpp
  ipc/embeddedSurfaceIPCPotentialEnergy.cpp
  ipc/broadPhase/spatialHashGrid.cpp
  ipc/broadPhase/surfaceIPCBroadPhase.cpp
  ipc/broadPhase/surfaceIPCSelfBroadPhase.cpp
  ipc/broadPhase/surfaceIPCExternalBroadPhase.cpp
  ipc/core/surfaceIPCSelfBarrierAssembler.cpp
  ipc/core/surfaceIPCExternalBarrierAssembler.cpp
  ipc/core/surfaceIPCBarrierKernels.cpp
  ipc/core/surfaceIPCMaxStep.cpp
  ipc/core/surfaceIPCCore.cpp
  ipc/geometry/ipcBarrier.cpp
  ipc/geometry/ipcCCD.cpp
  ipc/external/obstacleSurface.cpp
  ipc/external/obstaclePoseCache.cpp
  ipc/geometry/ipcDistancePrimitives.cpp
  ipc/geometry/ipcHessianProjection.cpp
  ipc/topology/surfaceIPCTopology.cpp
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
