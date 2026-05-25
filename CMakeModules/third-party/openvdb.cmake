if(TARGET openvdb_static OR TARGET OpenVDB::openvdb)
  return()
endif()

message(STATUS "Loading OpenVDB...")

pgo_dep_option(OPENVDB_BUILD_MAYA_PLUGIN BOOL OFF "Build OpenVDB Maya plugin")
pgo_dep_option(OPENVDB_ENABLE_UNINSTALL BOOL OFF "Adds a CMake uninstall target.")
pgo_dep_option(USE_HOUDINI BOOL OFF "Houdini")
pgo_dep_option(USE_MAYA BOOL OFF "Maya")
pgo_dep_option(USE_BLOSC BOOL OFF "Use Blosc compression")

if(WIN32)
  pgo_dep_option(USE_ZLIB BOOL OFF "Use ZLIB")
endif()

include(FetchContent)
FetchContent_Declare(
    openvdb
    URL https://github.com/AcademySoftwareFoundation/openvdb/archive/refs/tags/v12.1.1.zip
    EXCLUDE_FROM_ALL
    DOWNLOAD_EXTRACT_TIMESTAMP ON
)
pgo_fetch_make_available(openvdb)

if(DEFINED boost_SOURCE_DIR)
  file(GLOB PGO_OPENVDB_BOOST_INCLUDE_DIRS "${boost_SOURCE_DIR}/libs/*/include")
  foreach(tgt openvdb_static openvdb_shared)
    if(TARGET ${tgt})
      target_include_directories(${tgt} PRIVATE ${PGO_OPENVDB_BOOST_INCLUDE_DIRS})
    endif()
  endforeach()
endif()

message(STATUS "Done.")
