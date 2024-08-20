set(INTERPOLATIONCOORDINATES_HEADERS interpolationCoordinatesBase.h meanValueCoordinates.h GreenCoordinates.h barycentricCoordinates.h)
set(INTERPOLATIONCOORDINATES_SOURCES meanValueCoordinates.cpp GreenCoordinates.cpp barycentricCoordinates.cpp)
set(INTERPOLATIONCOORDINATES_DEPS volumetricMesh)

foreach(tgt ${INTERPOLATIONCOORDINATES_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "interpolationCoordinates is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(interpolationCoordinates "${INTERPOLATIONCOORDINATES_SOURCES}" "${INTERPOLATIONCOORDINATES_HEADERS}")

target_link_libraries(interpolationCoordinates PUBLIC ${INTERPOLATIONCOORDINATES_DEPS})