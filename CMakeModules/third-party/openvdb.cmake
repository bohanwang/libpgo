if(TARGET openvdb_static OR TARGET OpenVDB::openvdb)
  return()
endif()

message(STATUS "Loading OpenVDB...")

set(OPENVDB_BUILD_MAYA_PLUGIN OFF CACHE BOOL "Build OpenVDB Maya plugin" FORCE)
set(OPENVDB_ENABLE_UNINSTALL OFF CACHE BOOL "Adds a CMake uninstall target." FORCE)
set(USE_HOUDINI OFF CACHE BOOL "Houdini" FORCE)
set(USE_MAYA OFF CACHE BOOL "Maya" FORCE)

if(WIN32)
  set(USE_BLOSC CACHE BOOL "Maya" OFF)
  set(USE_ZLIB CACHE BOOL "Maya" OFF)
endif()

include(FetchContent)
FetchContent_Declare(
    openvdb
    URL https://github.com/AcademySoftwareFoundation/openvdb/archive/refs/tags/v12.1.1.zip
    EXCLUDE_FROM_ALL
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    FIND_PACKAGE_ARGS NAMES OpenVDB
)
FetchContent_MakeAvailable(openvdb)

message(STATUS "Done.")
