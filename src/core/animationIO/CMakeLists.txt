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

set(ANIMATIONIO_PUBLIC_DEPS pgoLogging mesh volumetricMesh configFileJSON interpolationCoordinates)
set(ANIMATIONIO_PRIVATE_DEPS Alembic::Alembic tiny_obj_loader)

foreach (tgt ${ANIMATIONIO_PUBLIC_DEPS} ${ANIMATIONIO_PRIVATE_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "animationIO is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

if(TARGET openvdb_static OR TARGET OpenVDB::openvdb)
  list(APPEND ANIMATIONIO_HEADERS stressFieldVDBExporter.h)
  list(APPEND ANIMATIONIO_SOURCES stressFieldVDBExporter.cpp)

  if(TARGET openvdb_static)
    list(APPEND ANIMATIONIO_PRIVATE_DEPS openvdb_static)
  else()
    list(APPEND ANIMATIONIO_PRIVATE_DEPS OpenVDB::openvdb)
  endif()

  message(STATUS "animationIO: stressFieldVDBExporter enabled (OpenVDB found).")
else()
  message(STATUS "animationIO: stressFieldVDBExporter disabled (OpenVDB not found; enable PGO_ENABLE_OPENVDB).")
endif()

add_libpgo_lib(animationIO "${ANIMATIONIO_SOURCES}" "${ANIMATIONIO_HEADERS}")

target_link_libraries(animationIO PUBLIC ${ANIMATIONIO_PUBLIC_DEPS})
target_link_libraries(animationIO PRIVATE ${ANIMATIONIO_PRIVATE_DEPS})
