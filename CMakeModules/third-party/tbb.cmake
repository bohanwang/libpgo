if(TARGET TBB::tbb)
  return()
endif()

message(STATUS "Loading tbb...")

set(TBB_DIR)

message(STATUS "conda: $ENV{CONDA_PREFIX}")

if(PGO_CHECK_CONDA AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
    message(STATUS "TBB_CMAKE:${CANDIDATE_PATH}")

    if(EXISTS ${CANDIDATE_PATH})
      set(TBB_DIR ${CANDIDATE_PATH})
    endif()
  elseif(WIN32)
    set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/Library/lib/cmake/TBB")

    if(EXISTS ${CANDIDATE_PATH})
      set(TBB_DIR ${CANDIDATE_PATH})
    endif()
  endif()
else()
  message(STATUS "TBB Not using conda")
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(EXISTS /opt/intel/oneapi/tbb/latest/lib/cmake/tbb)
      set(TBB_DIR /opt/intel/oneapi/tbb/latest/lib/cmake/tbb)
    endif()
  elseif(WIN32)
    if(EXISTS "C:/Program Files (x86)/Intel/oneAPI/tbb/latest/lib/cmake/tbb")
      set(TBB_DIR "C:/Program Files (x86)/Intel/oneAPI/tbb/latest/lib/cmake/tbb")
    endif()
  endif()
endif()

if(NOT TBB_DIR)
  include(FetchContent)
  FetchContent_Declare(
    tbb
    URL https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.12.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  FetchContent_GetProperties(tbb)

  if(NOT tbb_POPULATED)
    FetchContent_Populate(tbb)
  endif()

  set(TBB_TEST OFF CACHE BOOL "Enable testing" FORCE)
  set(TBB_EXAMPLES OFF CACHE BOOL "Enable examples")
  set(TBB_STRICT ON CACHE BOOL "Treat compiler warnings as errors" FORCE)
  set(TBB_WINDOWS_DRIVER OFF CACHE BOOL "Build as Universal Windows Driver (UWD)" FORCE)
  set(TBB_NO_APPCONTAINER OFF CACHE BOOL "Apply /APPCONTAINER:NO (for testing binaries for Windows Store)" FORCE)
  set(TBB4PY_BUILD OFF CACHE BOOL "Enable tbb4py build" FORCE)
  set(TBB_BUILD ON CACHE BOOL "Enable tbb build" FORCE)
  set(TBBMALLOC_BUILD ON CACHE BOOL "Enable tbbmalloc build" FORCE)
  set(TBB_CPF OFF CACHE BOOL "Enable preview features of the library" FORCE)
  set(TBB_FIND_PACKAGE OFF CACHE BOOL "Enable search for external oneTBB using find_package instead of build from sources" FORCE)
  set(TBB_DISABLE_HWLOC_AUTOMATIC_SEARCH OFF CACHE BOOL "Disable HWLOC automatic search by pkg-config tool" FORCE)
  set(TBB_ENABLE_IPO ON CACHE BOOL "Enable Interprocedural Optimization (IPO) during the compilation" FORCE)
  set(TBB_FUZZ_TESTING OFF CACHE BOOL "Enable fuzz testing" FORCE)
  set(TBB_INSTALL ON CACHE BOOL "Enable installation" FORCE)

  add_subdirectory(${tbb_SOURCE_DIR} ${tbb_BINARY_DIR})
  set(TBB_FOUND ON)
  mark_as_advanced(TBB_FOUND)
else()
  find_package(TBB CONFIG REQUIRED)
endif()

if(TARGET TBB::tbb)
  # fix the problem that some machine does not have release
  get_property(TBB_LIB TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE)

  if(TBB_LIB)
    message(STATUS "tbb lib: ${TBB_LIB}")
  else()
    get_property(TBB_LIB TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELWITHDEBINFO)
    if(TBB_LIB)
      set_target_properties(TBB::tbb PROPERTIES IMPORTED_LOCATION_RELEASE ${TBB_LIB})
    endif()

    get_property(TBB_LIB TARGET TBB::tbb PROPERTY IMPORTED_LOCATION_RELEASE)
    message(STATUS "tbb lib: ${TBB_LIB}")
  endif()
endif()

message(STATUS "Done.")
