if(NOT TARGET animationIO)
  message(STATUS "dumpStressVDB tool is not included. Missing: animationIO")
  return()
endif()

if(NOT (TARGET openvdb_static OR TARGET OpenVDB::openvdb))
  message(STATUS "dumpStressVDB tool is not included. Missing: OpenVDB (enable PGO_ENABLE_OPENVDB).")
  return()
endif()

add_libpgo_tools(dumpStressVDB dumpStressVDB.cpp "animationIO")
target_link_libraries(dumpStressVDB PRIVATE argparse::argparse)
