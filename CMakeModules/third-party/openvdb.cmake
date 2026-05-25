if(TARGET openvdb_static OR TARGET OpenVDB::openvdb)
  return()
endif()

message(STATUS "Loading OpenVDB...")

if(PGO_CHECK_CONDA AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
  if(WIN32)
    set(CANDIDATE_PREFIX "$ENV{CONDA_PREFIX}/Library")
  else()
    set(CANDIDATE_PREFIX "$ENV{CONDA_PREFIX}")
  endif()

  if(EXISTS "${CANDIDATE_PREFIX}")
    list(PREPEND CMAKE_PREFIX_PATH "${CANDIDATE_PREFIX}")
  endif()

  set(CANDIDATE_CONFIG_DIR "${CANDIDATE_PREFIX}/lib/cmake/OpenVDB")
  if(EXISTS "${CANDIDATE_CONFIG_DIR}/OpenVDBConfig.cmake")
    set(OpenVDB_DIR "${CANDIDATE_CONFIG_DIR}" CACHE PATH "OpenVDB CMake package directory" FORCE)
  elseif(EXISTS "${CANDIDATE_CONFIG_DIR}/FindOpenVDB.cmake")
    list(PREPEND CMAKE_MODULE_PATH "${CANDIDATE_CONFIG_DIR}")
  endif()
endif()

find_package(OpenVDB REQUIRED)

message(STATUS "Done.")
