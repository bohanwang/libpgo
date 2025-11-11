if(TARGET Alembic::Alembic)
  return()
endif()

message(STATUS "Loading alembic...")

set(USE_ARNOLD OFF CACHE BOOL "Include Arnold stuff" FORCE)
set(USE_BINARIES ON CACHE BOOL "Include binaries" FORCE)
set(USE_EXAMPLES OFF CACHE BOOL "Include examples" FORCE)
set(USE_HDF5 OFF CACHE BOOL "Include HDF5 stuff" FORCE)
set(USE_MAYA OFF CACHE BOOL "Include Maya stuff" FORCE)
set(USE_PRMAN OFF CACHE BOOL "Include PRMan stuff" FORCE)
set(USE_PYALEMBIC OFF CACHE BOOL "Include PyAlembic stuff" FORCE)
set(USE_STATIC_BOOST OFF CACHE BOOL "Build with static Boost libs" FORCE)
set(USE_STATIC_HDF5 OFF CACHE BOOL "Build with static HDF5 libs" FORCE)
set(USE_TESTS OFF CACHE BOOL "Include Alembic tests" FORCE)
set(ALEMBIC_BUILD_LIBS ON CACHE BOOL "Build library, if off use external alembic libs" FORCE)
set(ALEMBIC_SHARED_LIBS ON CACHE BOOL "Build shared libraries" FORCE)
set(ALEMBIC_DEBUG_WARNINGS_AS_ERRORS ON CACHE BOOL "In debug mode build with warnings as errors" FORCE)

include(FetchContent)
FetchContent_Declare(
  alembic
  URL https://github.com/alembic/alembic/archive/refs/tags/1.8.9.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES alembic
)

FetchContent_MakeAvailable(alembic)

message(STATUS "Done.")





