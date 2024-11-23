if(TARGET geogram::geogram)
  return()
endif()

message(STATUS "Loading geogram...")

set(GEOGRAM_SUB_BUILD ON CACHE BOOL "" FORCE)
set(GEOGRAM_WITH_HLBFGS ON CACHE BOOL "Non-linear solver (Yang Liu's HLBFGS)" FORCE)

include(FetchContent)
FetchContent_Declare(
  geogram
  URL https://github.com/BrunoLevy/geogram/releases/download/v1.9.0/geogram_1.9.0.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES geogram
)

FetchContent_GetProperties(geogram)
if(NOT geogram_POPULATED)
  FetchContent_Populate(geogram)
endif()

set(MODIFIED_FILE "${CMAKE_SOURCE_DIR}/CMakeModules/patches/geogram.cmake")
set(TARGET_FILE "${geogram_SOURCE_DIR}/CMakeLists.txt")

# Read in the content
file(READ "${MODIFIED_FILE}" content)

# Write the modified content back to the file
file(WRITE "${TARGET_FILE}" "${content}")

add_subdirectory(${geogram_SOURCE_DIR} ${geogram_BINARY_DIR} EXCLUDE_FROM_ALL)

set_target_properties(geogram PROPERTIES CXX_STANDARD 14)

message(STATUS "Done.")