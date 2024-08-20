set(GEOGRAM_INTERFACE_HEADERS geogramInterface.h)
set(GEOGRAM_INTERFACE_SOURCES geogramInterface.cpp)
set(GEOGRAM_INTERFACE_DEPS pgoLogging eigenSupport mesh geogram)

foreach (tgt ${GEOGRAM_INTERFACE_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "geogramInterface is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(geogramInterface "${GEOGRAM_INTERFACE_SOURCES}" "${GEOGRAM_INTERFACE_HEADERS}")

target_link_libraries(geogramInterface PUBLIC ${GEOGRAM_INTERFACE_DEPS})