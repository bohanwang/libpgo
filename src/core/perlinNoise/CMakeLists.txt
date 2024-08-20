set(PERLINNOISE_HEADERS
  perlinNoise.h
  perlinWind.h
)

set(PERLINNOISE_SOURCES
  perlinNoise.cpp
  perlinWind.cpp
)

set(PERLINNOISE_DEPS eigenSupport)

foreach(tgt ${PERLINNOISE_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "perlinNoise is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(perlinNoise "${PERLINNOISE_SOURCES}" "${PERLINNOISE_HEADERS}")

target_link_libraries(perlinNoise PUBLIC ${PERLINNOISE_DEPS})