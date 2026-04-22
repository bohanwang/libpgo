set(PROFILING_HEADERS
  scopedProfileSection.h
)

set(PROFILING_SOURCES
  scopedProfileSection.cpp
)

add_libpgo_lib(profiling "${PROFILING_SOURCES}" "${PROFILING_HEADERS}")
