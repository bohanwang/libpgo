if(TARGET Eigen3::Eigen)
else()
  message(STATUS "Loading eigen...")

  set(BUILD_TESTING OFF CACHE BOOL "eigen build test" FORCE)

  include(FetchContent)
  FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    EXCLUDE_FROM_ALL
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    FIND_PACKAGE_ARGS NAMES Eigen3
  )

  FetchContent_MakeAvailable(eigen)

  message(STATUS "Done.")
endif()

if(TARGET MKL::MKL)
  get_property(aliased_target TARGET Eigen3::Eigen PROPERTY ALIASED_TARGET)
  if("${aliased_target}" STREQUAL "")
    set(REAL_TGT Eigen3::Eigen)
  else()
    set(REAL_TGT ${aliased_target})
  endif()
  
  target_link_libraries(${REAL_TGT} INTERFACE MKL::MKL)
  target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_DONT_PARALLELIZE)
  target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_USE_MKL_ALL)
endif()
