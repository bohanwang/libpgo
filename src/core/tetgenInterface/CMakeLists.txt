set(TETGEN_INTERFACE_HEADERS tetgenInterface.h)
set(TETGEN_INTERFACE_SOURCES tetgenInterface.cpp)
set(TETGEN_INTERFACE_DEPS pgoLogging eigenSupport mesh tetgen)

foreach (tgt ${TETGEN_INTERFACE_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "tetgenInterface is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(tetgenInterface "${TETGEN_INTERFACE_SOURCES}" "${TETGEN_INTERFACE_HEADERS}")

target_link_libraries(tetgenInterface PUBLIC ${TETGEN_INTERFACE_DEPS})