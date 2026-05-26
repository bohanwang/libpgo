set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

set(RUN_IPC_SIM_CORE_SOURCES
  cli/cliLogging.cpp
  io/volumeMeshIO.cpp
  setup/femSetup.cpp
  setup/setup.cpp
  setup/setupCommon.cpp
  setup/attachmentSetup.cpp
  setup/floorSetup.cpp
  setup/surfacePressureSetup.cpp
  setup/obstacleSetup.cpp
  setup/shellSetup.cpp
  setup/volumeSetup.cpp
  setup/legacySetup.cpp
  app/app.cpp
  app/config.cpp
  app/output.cpp
  app/logging.cpp
  app/session.cpp
  app/loop.cpp
  solver/solverRunner.cpp
  solver/staticSolve.cpp
  contact/ipcContactBackend.cpp
  contact/legacyPenaltyContact.cpp
)

add_library(runIPCSimCore STATIC ${RUN_IPC_SIM_CORE_SOURCES})
target_link_libraries(runIPCSimCore PUBLIC ${AVAILABLE_LIBS})
target_link_libraries(runIPCSimCore PRIVATE compilation_flag)
target_include_directories(runIPCSimCore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET runIPCSimCore PROPERTY FOLDER "tools/internal")

if(PGO_RELEASE_MODE_DEBUG)
  target_link_libraries(runIPCSimCore PRIVATE compilation_flag_for_debug)
endif()

add_libpgo_tools(runIPCSim "runIPCSim.cpp;cli/cli.cpp" "")
target_link_libraries(runIPCSim PRIVATE runIPCSimCore argparse::argparse)
