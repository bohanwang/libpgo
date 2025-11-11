set(ANIMATIONIO_HEADERS
  animationLoader.h
  abcWriter.h
)

set(ANIMATIONIO_SOURCES
  animationLoader.cpp
  abcWriter.cpp
)

if(NOT TARGET Alembic::Alembic)
  message(STATUS "animationIO is not included. Missing: Alembic::Alembic")
  return()
endif()

set(ANIMATIONIO_DEPS pgoLogging mesh volumetricMesh configFileJSON interpolationCoordinates Alembic::Alembic)

foreach (tgt ${ANIMATIONIO_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "animationIO is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(animationIO "${ANIMATIONIO_SOURCES}" "${ANIMATIONIO_HEADERS}")

target_link_libraries(animationIO PUBLIC ${ANIMATIONIO_DEPS})
