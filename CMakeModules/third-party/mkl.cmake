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

  # oneMKL 2025.x checks this variable even when the imported TBB target
  # already exists.
  if(TARGET TBB::tbb)
    set(TBB_tbb_FOUND TRUE)
  endif()

  find_package(MKL CONFIG REQUIRED)
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
      set(
        _pgo_mkl_dispatch_library
        "${MKL_ROOT}/lib/libmkl_${_pgo_mkl_dispatch_name}.so.2"
      )
      if(NOT EXISTS "${_pgo_mkl_dispatch_library}")
        message(
          FATAL_ERROR
          "Required oneMKL dispatch library is missing: "
          "${_pgo_mkl_dispatch_library}"
        )
      endif()
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
