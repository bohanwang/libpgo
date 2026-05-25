if(TARGET Eigen3::Eigen)
else()
  message(STATUS "Loading eigen...")

  pgo_dep_option(BUILD_TESTING BOOL OFF "eigen build test")
  pgo_dep_option(BUILD_EXAMPLES BOOL OFF "eigen build examples")
  pgo_dep_option(EIGEN_BUILD_DOC BOOL OFF "eigen build doc")
  pgo_dep_option(EIGEN_BUILD_CMAKE_PACKAGE BOOL ON "eigen build cmake package")

  include(FetchContent)
  FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    # GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    # GIT_TAG fb2fca90be39783f76ba05b521360afeeda265f2
    EXCLUDE_FROM_ALL
    DOWNLOAD_EXTRACT_TIMESTAMP ON
  )

  pgo_fetch_make_available(Eigen3)

  message(STATUS "Done.")
endif()

get_property(aliased_target TARGET Eigen3::Eigen PROPERTY ALIASED_TARGET)
if("${aliased_target}" STREQUAL "")
  set(REAL_TGT Eigen3::Eigen)
else()
  set(REAL_TGT ${aliased_target})
endif()

if(TARGET MKL::MKL)
  target_link_libraries(${REAL_TGT} INTERFACE MKL::MKL)
  target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_DONT_PARALLELIZE)
  target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_USE_MKL_ALL)
  target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_MKL_NO_DIRECT_CALL)
endif()

target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_MAX_ALIGN_BYTES=32)
