if(TARGET pgoLogging)
  add_libpgo_tools(computeStressFieldStats computeStressFieldStats.cpp "pgoLogging")
  target_link_libraries(computeStressFieldStats PRIVATE argparse::argparse nlohmann_json::nlohmann_json fmt::fmt)
endif()

if(NOT TARGET implicitSurface)
  message(STATUS "generateFBMSUnionSurface tool is not included. Missing: implicitSurface")
  return()
endif()

set(GENERATE_FBMS_UNION_SURFACE_ADAPTER_SOURCES
  generateFBMSUnionSurface.cpp
  generateFBMSUnionSurfaceAdapter.cpp
  generateFBMSUnionSurfaceIO.cpp
  generateFBMSUnionSurfaceOptions.cpp)

set(GENERATE_FBMS_UNION_SURFACE_ADAPTER_HEADERS
  generateFBMSUnionSurfaceAdapter.h
  generateFBMSUnionSurfaceIO.h
  generateFBMSUnionSurfaceOptions.h
  generateFBMSUnionSurfaceVolume.h)

add_libpgo_tools(generateFBMSUnionSurface
  "${GENERATE_FBMS_UNION_SURFACE_ADAPTER_SOURCES}"
  "")

target_sources(generateFBMSUnionSurface PRIVATE
  ${GENERATE_FBMS_UNION_SURFACE_ADAPTER_HEADERS})

target_link_libraries(generateFBMSUnionSurface PRIVATE
  implicitSurface
  argparse::argparse)
