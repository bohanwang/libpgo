if(TARGET TBB::tbb)
  return()
endif()

message(STATUS "Loading tbb...")

find_package(TBB CONFIG REQUIRED)
message(STATUS "Using external TBB package: ${TBB_DIR}")

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
