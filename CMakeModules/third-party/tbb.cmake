if(TARGET TBB::tbb)
  return()
endif()

message(STATUS "Loading tbb...")

if(PGO_CHECK_CONDA AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
  if(WIN32)
    set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/Library/lib/cmake/TBB")
  else()
    set(CANDIDATE_PATH "$ENV{CONDA_PREFIX}/lib/cmake/TBB")
  endif()

  if(EXISTS "${CANDIDATE_PATH}/TBBConfig.cmake")
    set(TBB_DIR "${CANDIDATE_PATH}")
  endif()
endif()

find_package(TBB CONFIG REQUIRED)

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
