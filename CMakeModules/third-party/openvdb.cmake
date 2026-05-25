if(TARGET openvdb_static OR TARGET OpenVDB::openvdb)
  return()
endif()

message(STATUS "Loading OpenVDB...")

if(PGO_CHECK_CONDA AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
  if(WIN32)
    set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/Library/lib/cmake/OpenVDB")
  else()
    set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/lib/cmake/OpenVDB")
  endif()

  if(EXISTS "${CANDIDATE_PATH}/FindOpenVDB.cmake")
    list(PREPEND CMAKE_MODULE_PATH "${CANDIDATE_PATH}")
  endif()
endif()

find_package(OpenVDB REQUIRED)

message(STATUS "Done.")
