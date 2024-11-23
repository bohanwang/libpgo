if(TARGET MKL::MKL)
  return()
endif()

if(PGO_USE_MKL)
  message(STATUS "Searching for MKL")
  
  set(MKL_THREADING tbb_thread)
  set(MKL_INTERFACE lp64)

  if(PGO_MKL_LINK_DYNAMIC)
    set(MKL_LINK dynamic)
  else()
    set(MKL_LINK static)
  endif()
  
  find_package(MKL CONFIG)

  if(TARGET MKL::MKL)
    set(HAS_MKL 1)
  else()
    set(HAS_MKL 0)
    message(FATAL_ERROR "Cannot find MKL")
  endif()
else()
  message(STATUS "Not using MKL")
  set(HAS_MKL 0)
endif()