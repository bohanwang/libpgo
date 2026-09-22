if(TARGET MKL::MKL)
  return()
endif()

if(PGO_USE_MKL)
  message(STATUS "Searching for MKL")

  # Project OpenMP is independent of oneMKL's internal threading backend.
  # Keep Windows and Linux on oneTBB and reject cache overrides that could
  # silently pull in an Intel/GNU OpenMP threading runtime.
  set(MKL_THREADING tbb_thread CACHE STRING "oneMKL threading layer" FORCE)
  set(MKL_INTERFACE lp64 CACHE STRING "oneMKL integer interface" FORCE)

  if(PGO_MKL_LINK_DYNAMIC)
    set(MKL_LINK dynamic)
  else()
    set(MKL_LINK static)
  endif()

  # oneMKL 2025.x checks this variable even when the imported TBB target
  # already exists.
  if(TARGET TBB::tbb)
    set(TBB_tbb_FOUND TRUE)
  endif()

  find_package(MKL CONFIG REQUIRED)
  target_link_options(MKL::MKL INTERFACE "LINKER:--no-as-needed")
  message(STATUS "Using external MKL package: ${MKL_DIR}")

  # oneMKL loads its CPU dispatch kernels with dlopen(), so repair tools cannot
  # discover them from the normal ELF dependency closure. Keep the complete
  # Linux dispatch set as explicit NEEDED entries so auditwheel can vendor
  # them. The wheel workflow restores only these filenames after repair because
  # oneMKL's embedded dlopen() strings cannot follow auditwheel's hash renames.
  if(PGO_MKL_LINK_DYNAMIC AND UNIX AND NOT APPLE)
    set(_pgo_mkl_dispatch_libraries)
    foreach(_pgo_mkl_dispatch_name
        def
        mc3
        avx2
        avx512
        vml_def
        vml_cmpt
        vml_mc3
        vml_avx2
        vml_avx512)
      # The SONAME version changes between oneMKL releases (.so.2 up to 2025,
      # .so.3 from 2026), so discover the versioned file instead of assuming it.
      file(GLOB _pgo_mkl_dispatch_candidates
        "${MKL_ROOT}/lib/libmkl_${_pgo_mkl_dispatch_name}.so.[0-9]*")
      list(LENGTH _pgo_mkl_dispatch_candidates _pgo_mkl_dispatch_count)
      if(_pgo_mkl_dispatch_count EQUAL 0)
        message(
          FATAL_ERROR
          "Required oneMKL dispatch library is missing: "
          "${MKL_ROOT}/lib/libmkl_${_pgo_mkl_dispatch_name}.so.<version>"
        )
      endif()
      list(SORT _pgo_mkl_dispatch_candidates COMPARE NATURAL ORDER DESCENDING)
      list(GET _pgo_mkl_dispatch_candidates 0 _pgo_mkl_dispatch_library)
      list(APPEND _pgo_mkl_dispatch_libraries "${_pgo_mkl_dispatch_library}")
    endforeach()
    target_link_libraries(
      MKL::MKL
      INTERFACE
        "-Wl,--no-as-needed"
        ${_pgo_mkl_dispatch_libraries}
        "-Wl,--as-needed"
    )
  endif()

  set(HAS_MKL 1)
else()
  message(STATUS "Not using MKL")
  set(HAS_MKL 0)
endif()
