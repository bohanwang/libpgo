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

  # OpenVDB's conda package is built against conda Boost. Keep its dependency
  # checks inside the same conda prefix.
  set(BOOST_ROOT "${CANDIDATE_PREFIX}" CACHE PATH "Boost prefix for OpenVDB" FORCE)
  set(BOOST_INCLUDEDIR "${CANDIDATE_PREFIX}/include" CACHE PATH "Boost include directory for OpenVDB" FORCE)
  set(Boost_INCLUDE_DIR "${CANDIDATE_PREFIX}/include" CACHE PATH "Boost include directory for OpenVDB" FORCE)
  set(Boost_INCLUDE_DIRS "${CANDIDATE_PREFIX}/include" CACHE STRING "Boost include directories for OpenVDB" FORCE)
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_STATIC_LIBS OFF CACHE BOOL "Use conda shared Boost libraries for OpenVDB" FORCE)
  set(OPENVDB_USE_STATIC_LIBS OFF)
  set(OPENVDB_USE_STATIC_LIBS OFF CACHE BOOL "Use conda shared OpenVDB libraries" FORCE)
  set(Boost_NO_SYSTEM_PATHS ON CACHE BOOL "Restrict OpenVDB Boost lookup to conda" FORCE)
  set(Boost_NO_BOOST_CMAKE OFF CACHE BOOL "Allow OpenVDB to find conda Boost config" FORCE)
  set(BOOST_LIBRARYDIR "${CANDIDATE_PREFIX}/lib" CACHE PATH "Boost library directory for OpenVDB" FORCE)

  set(CANDIDATE_CONFIG_DIR "${CANDIDATE_PREFIX}/lib/cmake/OpenVDB")
  if(EXISTS "${CANDIDATE_CONFIG_DIR}/OpenVDBConfig.cmake")
    set(OpenVDB_DIR "${CANDIDATE_CONFIG_DIR}" CACHE PATH "OpenVDB CMake package directory" FORCE)
  elseif(EXISTS "${CANDIDATE_CONFIG_DIR}/FindOpenVDB.cmake")
    list(PREPEND CMAKE_MODULE_PATH "${CANDIDATE_CONFIG_DIR}")
  endif()
endif()

find_package(OpenVDB REQUIRED)

message(STATUS "Done.")
