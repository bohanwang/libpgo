set(AVAILABLE_LIBS mesh)
set(MESH_QUALITY_CHECK_HAS_CGAL 0)

if(TARGET cgalInterface)
  list(APPEND AVAILABLE_LIBS cgalInterface)
  set(MESH_QUALITY_CHECK_HAS_CGAL 1)
endif()

if(NOT TARGET mesh)
  message(STATUS "meshQualityCheck tool is not included. Missing: mesh")
  return()
endif()

add_libpgo_tools(meshQualityCheck "meshQualityCheck.cpp;surfaceQuality.cpp" "${AVAILABLE_LIBS}")
target_compile_definitions(meshQualityCheck PRIVATE MESH_QUALITY_CHECK_HAS_CGAL=${MESH_QUALITY_CHECK_HAS_CGAL})
target_link_libraries(meshQualityCheck PRIVATE argparse::argparse nlohmann_json::nlohmann_json)
