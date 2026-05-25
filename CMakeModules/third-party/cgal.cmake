if(TARGET CGAL::CGAL)
  return()
endif()

message(STATUS "Downloading CGAL...")
include(FetchContent)

FetchContent_Declare(
  CGAL
  URL https://github.com/CGAL/cgal/releases/download/v6.0.1/CGAL-6.0.1-library.zip
  URL_HASH SHA256=7b0bb231d57261491722b7a0950f8026e17a08ae4a93315495cfb9d91faa31e3
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_populate_compat(CGAL "CGAL source tree is patched before find_package(CGAL CONFIG ... NO_DEFAULT_PATH)")

set(CGAL_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/cgal-src")
set(CGAL_BINARY_DIR "${CMAKE_BINARY_DIR}/_deps/cgal-build")

set(MODIFIED_FILE "${CMAKE_SOURCE_DIR}/CMakeModules/patches/CGAL_SetupBoost.cmake")
set(CGAL_SETUP_BOOST_FILE "${CGAL_SOURCE_DIR}/cmake/modules/CGAL_SetupBoost.cmake")
if(NOT EXISTS "${CGAL_SETUP_BOOST_FILE}")
  message(FATAL_ERROR "CGAL setup file not found: ${CGAL_SETUP_BOOST_FILE}")
endif()

file(READ "${MODIFIED_FILE}" content)
file(WRITE "${CGAL_SETUP_BOOST_FILE}" "${content}")

# gmp for windows
if(WIN32)
  # message(STATUS "Downloading gmp...")
  # set(GMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/cgal-gmp.zip")
  # file(DOWNLOAD https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1-win64-auxiliary-libraries-gmp-mpfr.zip ${GMP_FILE})
  # file(ARCHIVE_EXTRACT INPUT ${GMP_FILE} DESTINATION ${cgal_SOURCE_DIR})
  if(PGO_CHECK_CONDA AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
    set(TGT_FILE "$ENV{CONDA_PREFIX}/Library/bin/gmp-10.dll")

    if(NOT EXISTS ${TGT_FILE})
      file(COPY_FILE "${CMAKE_SOURCE_DIR}/third-party/gmp-msvc/release/gmp-10.dll" "${TGT_FILE}")
    endif()

    set(TGT_FILE "$ENV{CONDA_PREFIX}/Library/bin/gmpxx-4.dll")

    if(NOT EXISTS ${TGT_FILE})
      file(COPY_FILE "${CMAKE_SOURCE_DIR}/third-party/gmp-msvc/release/gmpxx-4.dll" "${TGT_FILE}")
    endif()

    set(TGT_FILE "$ENV{CONDA_PREFIX}/Library/bin/mpfr-6.dll")

    if(NOT EXISTS ${TGT_FILE})
      file(COPY_FILE "${CMAKE_SOURCE_DIR}/third-party/mpfr-msvc/release/mpfr-6.dll" "${TGT_FILE}")
    endif()
  endif()
endif()

# CGAL is header-only here but hard-depends on the native GMP/MPFR (+GMPXX)
# libraries. On macOS the Homebrew prefix is added to the global search paths in
# the top-level CMakeLists, so CGAL's own find_package(GMP) can resolve them.
#
# Make CGAL's FindGMP/FindMPFR/FindGMPXX modules available, then locate the
# libraries so they are found before (and reused by) find_package(CGAL). These are
# best-effort (not REQUIRED): on some platforms GMP is only discoverable through
# CGAL's own auxiliary-dir hint, which is set during find_package(CGAL). CGAL's
# internal CGAL_SetupGMP does the authoritative REQUIRED find, so a miss here is
# harmless.
list(APPEND CMAKE_MODULE_PATH "${CGAL_SOURCE_DIR}/cmake/modules")
find_package(GMP QUIET)
find_package(MPFR QUIET)
find_package(GMPXX QUIET)

pgo_dep_option(CGAL_WITH_GMPXX BOOL ON "Enable CGAL GMPXX support")
pgo_dep_option(CGAL_ENABLE_TESTING BOOL OFF "disable testing")
find_package(CGAL CONFIG COMPONENTS Core REQUIRED PATHS ${CGAL_SOURCE_DIR} NO_DEFAULT_PATH)

message(STATUS "cgal module path: ${CGAL_MODULES_DIR}")
message(STATUS "TBB: ${TBB_FOUND}")
include("${CGAL_MODULES_DIR}/CGAL_TBB_support.cmake")
