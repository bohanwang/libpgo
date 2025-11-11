if(TARGET Eigen3::Eigen)
else()
  message(STATUS "Loading eigen...")

  set(BUILD_TESTING OFF CACHE BOOL "eigen build test" FORCE)
  set(EIGEN_BUILD_CMAKE_PACKAGE ON CACHE BOOL "eigen build cmake package" FORCE)

  include(FetchContent)
  FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    # GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    # GIT_TAG fb2fca90be39783f76ba05b521360afeeda265f2
    EXCLUDE_FROM_ALL
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    FIND_PACKAGE_ARGS NAMES Eigen3
  )

  FetchContent_MakeAvailable(Eigen3)

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

get_target_property(flags compilation_flag INTERFACE_COMPILE_OPTIONS)
message(STATUS "Eigen3 compilation flags: ${flags}")
target_compile_options(${REAL_TGT} INTERFACE ${flags})
target_compile_definitions(${REAL_TGT} INTERFACE EIGEN_MAX_ALIGN_BYTES=32)