if(TARGET TBB::tbb)
  return()
endif()

message(STATUS "Loading tbb...")

FetchContent_Declare(
  tbb
  URL https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.12.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

pgo_dep_option(TBB_TEST BOOL OFF "Enable testing")
pgo_dep_option(TBB_EXAMPLES BOOL OFF "Enable examples")
pgo_dep_option(TBB_STRICT BOOL ON "Treat compiler warnings as errors")
pgo_dep_option(TBB_WINDOWS_DRIVER BOOL OFF "Build as Universal Windows Driver (UWD)")
pgo_dep_option(TBB_NO_APPCONTAINER BOOL OFF "Apply /APPCONTAINER:NO (for testing binaries for Windows Store)")
pgo_dep_option(TBB4PY_BUILD BOOL OFF "Enable tbb4py build")
pgo_dep_option(TBB_BUILD BOOL ON "Enable tbb build")
pgo_dep_option(TBBMALLOC_BUILD BOOL ON "Enable tbbmalloc build")
pgo_dep_option(TBB_CPF BOOL OFF "Enable preview features of the library")
pgo_dep_option(TBB_FIND_PACKAGE BOOL OFF "Enable search for external oneTBB using find_package instead of build from sources")
pgo_dep_option(TBB_DISABLE_HWLOC_AUTOMATIC_SEARCH BOOL OFF "Disable HWLOC automatic search by pkg-config tool")
pgo_dep_option(TBB_ENABLE_IPO BOOL ON "Enable Interprocedural Optimization (IPO) during the compilation")
pgo_dep_option(TBB_FUZZ_TESTING BOOL OFF "Enable fuzz testing")
pgo_dep_option(TBB_INSTALL BOOL ON "Enable installation")

pgo_fetch_make_available(tbb)
set(TBB_FOUND ON)
mark_as_advanced(TBB_FOUND)

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
