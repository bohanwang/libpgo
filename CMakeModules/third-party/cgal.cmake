if(TARGET CGAL::CGAL)
  return()
endif()

message(STATUS "Downloading CGAL...")
include(FetchContent)
FetchContent_Declare(
  cgal
  URL https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1-library.tar.xz
  URL_HASH SHA256=6640f91dad9956764db27440932ac3a2f6269bf8066690f639dd7e995e3e1925
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES CGAL COMPONENTS Core
)

FetchContent_GetProperties(cgal)

if(NOT cgal_POPULATED)
  FetchContent_Populate(cgal)
endif()

set(MODIFIED_FILE "${CMAKE_SOURCE_DIR}/CMakeModules/patches/CGAL_SetupBoost.cmake")
set(TARGET_FILE "${cgal_SOURCE_DIR}/cmake/modules/CGAL_SetupBoost.cmake")

# Read in the content
file(READ "${MODIFIED_FILE}" content)

# Write the modified content back to the file
file(WRITE "${TARGET_FILE}" "${content}")

# gmp for windows
if(WIN32)
  # message(STATUS "Downloading gmp...")
  # set(GMP_FILE "${CMAKE_CURRENT_BINARY_DIR}/cgal-gmp.zip")
  # file(DOWNLOAD https://github.com/CGAL/cgal/releases/download/v5.6.1/CGAL-5.6.1-win64-auxiliary-libraries-gmp-mpfr.zip ${GMP_FILE})
  # file(ARCHIVE_EXTRACT INPUT ${GMP_FILE} DESTINATION ${cgal_SOURCE_DIR})
  if("$ENV{CONDA_PREFIX}" STREQUAL "")
  else()
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

set(CGAL_WITH_GMPXX ON CACHE BOOL "" FORCE)
find_package(CGAL CONFIG COMPONENTS Core REQUIRED PATHS ${cgal_SOURCE_DIR} NO_DEFAULT_PATH)

message(STATUS "cgal m: ${CGAL_MODULES_DIR}")
message(STATUS "TBB: ${TBB_FOUND}")

include("${CGAL_MODULES_DIR}/CGAL_TBB_support.cmake")

message(STATUS "Done.")

# FetchContent_GetProperties(cgal)

# message(STATUS ${cgal_SOURCE_DIR})

# set(CGAL_SOURCE_DIR ${cgal_SOURCE_DIR})

# FetchContent_MakeAvailable(cgal)

# message(STATUS "cgal m: ${CGAL_MODULES_DIR}")
# include("${CGAL_MODULES_DIR}/CGAL_TBB_support.cmake")

# message(STATUS "Done.")