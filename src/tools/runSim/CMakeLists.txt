set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

set(RUN_IPC_SIM_CORE_SOURCES
  runSimCliLogging.cpp
  runSimVolumeMeshIO.cpp
  runSimFEMSetup.cpp
  runIPCSimSetup.cpp
  runIPCSimSetupCommon.cpp
  runIPCSimAttachmentSetup.cpp
  runIPCSimFloorSetup.cpp
  runIPCSimSurfacePressureSetup.cpp
  runIPCSimObstacleSetup.cpp
  runIPCSimShellSetup.cpp
  runIPCSimVolumeSetup.cpp
  runIPCSimLegacySetup.cpp
  runIPCSimApp.cpp
  runIPCSimConfig.cpp
  runIPCSimOutput.cpp
  runIPCSimLogging.cpp
  runIPCSimSession.cpp
  runIPCSimLoop.cpp
  runIPCSimStaticSolve.cpp
  runIPCSimIpcContactBackend.cpp
  runIPCSimLegacyPenaltyContact.cpp
)

add_library(runIPCSimCore STATIC ${RUN_IPC_SIM_CORE_SOURCES})
target_link_libraries(runIPCSimCore PUBLIC ${AVAILABLE_LIBS})
target_link_libraries(runIPCSimCore PRIVATE compilation_flag)
target_include_directories(runIPCSimCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET runIPCSimCore PROPERTY FOLDER "tools/internal")

if(PGO_RELEASE_MODE_DEBUG)
  target_link_libraries(runIPCSimCore PRIVATE compilation_flag_for_debug)
endif()

add_libpgo_tools(runIPCSim "runIPCSim.cpp;runIPCSimCli.cpp" "")
target_link_libraries(runIPCSim PRIVATE runIPCSimCore argparse::argparse)
