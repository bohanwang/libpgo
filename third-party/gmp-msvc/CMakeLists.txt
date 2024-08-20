set(GMP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Include path to GMP" FORCE)
set(GMP_LIBRARIES debug "${CMAKE_CURRENT_SOURCE_DIR}/debug/gmp.lib" optimized "${CMAKE_CURRENT_SOURCE_DIR}/release/gmp.lib" CACHE PATH "Lib path to GMP" FORCE)

set(GMPXX_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Include path to GMPXX" FORCE)
set(GMPXX_LIBRARIES debug "${CMAKE_CURRENT_SOURCE_DIR}/debug/gmpxx.lib" optimized "${CMAKE_CURRENT_SOURCE_DIR}/release/gmpxx.lib" CACHE PATH "Lib path to GMPXX" FORCE)
