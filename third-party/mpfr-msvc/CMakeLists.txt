set(MPFR_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Include path to MPFR" FORCE)
set(MPFR_LIBRARIES debug "${CMAKE_CURRENT_SOURCE_DIR}/debug/mpfr.lib" optimized "${CMAKE_CURRENT_SOURCE_DIR}/release/mpfr.lib" CACHE PATH "Lib path to MPFR" FORCE)
