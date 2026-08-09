if(TARGET Eigen3::Eigen)
else()
  message(STATUS "Loading eigen...")

  set(BUILD_TESTING OFF CACHE BOOL "eigen build test" FORCE)
  set(BUILD_EXAMPLES OFF CACHE BOOL "eigen build examples" FORCE)
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "eigen build documentation" FORCE)
  set(EIGEN_BUILD_CMAKE_PACKAGE ON CACHE BOOL "eigen build cmake package" FORCE)

  include(FetchContent)
  FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
    # GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    # GIT_TAG fb2fca90be39783f76ba05b521360afeeda265f2
    OVERRIDE_FIND_PACKAGE
    EXCLUDE_FROM_ALL
    DOWNLOAD_EXTRACT_TIMESTAMP ON
  )

  # Eigen 3.4 configures its optional BLAS target even when BUILD_TESTING is
  # disabled. On Windows, CheckLanguage may pick up an unrelated MinGW
  # gfortran from PATH and then try to combine it with the active MSVC
  # toolchain. libpgo does not use Eigen's BLAS target, so prevent only that
  # optional probe while preserving an explicitly configured Fortran compiler.
  if(WIN32 AND NOT DEFINED CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER NOTFOUND)
    set(_pgo_suppressed_eigen_fortran_probe TRUE)
  endif()

  FetchContent_MakeAvailable(Eigen3)

  if(_pgo_suppressed_eigen_fortran_probe)
    unset(CMAKE_Fortran_COMPILER)
    unset(_pgo_suppressed_eigen_fortran_probe)
  endif()

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
