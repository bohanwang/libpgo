set(LIBIGL_INTERFACE_HEADERS libiglInterface.h)
set(LIBIGL_INTERFACE_SOURCES libiglInterface.cpp)
set(LIBIGL_INTERFACE_DEPS pgoLogging eigenSupport mesh igl::core igl_copyleft::cgal)

foreach (tgt ${LIBIGL_INTERFACE_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "libigl is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(libiglInterface "${LIBIGL_INTERFACE_SOURCES}" "${LIBIGL_INTERFACE_HEADERS}")

target_link_libraries(libiglInterface PUBLIC ${LIBIGL_INTERFACE_DEPS})