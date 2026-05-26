if(TARGET Alembic::Alembic)
  return()
endif()

message(STATUS "Loading alembic...")

pgo_dep_option(USE_ARNOLD BOOL OFF "Include Arnold stuff")
pgo_dep_option(USE_BINARIES BOOL OFF "Include binaries")
pgo_dep_option(USE_EXAMPLES BOOL OFF "Include examples")
pgo_dep_option(USE_HDF5 BOOL OFF "Include HDF5 stuff")
pgo_dep_option(USE_MAYA BOOL OFF "Include Maya stuff")
pgo_dep_option(USE_PRMAN BOOL OFF "Include PRMan stuff")
pgo_dep_option(USE_PYALEMBIC BOOL OFF "Include PyAlembic stuff")
pgo_dep_option(USE_STATIC_BOOST BOOL OFF "Build with static Boost libs")
pgo_dep_option(USE_STATIC_HDF5 BOOL OFF "Build with static HDF5 libs")
pgo_dep_option(USE_TESTS BOOL OFF "Include Alembic tests")
pgo_dep_option(ALEMBIC_BUILD_LIBS BOOL ON "Build library, if off use external alembic libs")
pgo_dep_option(ALEMBIC_SHARED_LIBS BOOL OFF "Build shared libraries")
pgo_dep_option(ALEMBIC_DEBUG_WARNINGS_AS_ERRORS BOOL ON "In debug mode build with warnings as errors")

include(FetchContent)
FetchContent_Declare(
  alembic
  URL https://github.com/alembic/alembic/archive/refs/tags/1.8.9.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(alembic)

message(STATUS "Done.")




