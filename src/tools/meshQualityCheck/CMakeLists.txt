if(NOT TARGET mesh)
  message(STATUS "meshQualityCheck tool is not included. Missing: mesh")
  return()
endif()

add_libpgo_tools(meshQualityCheck "meshQualityCheck.cpp;surfaceQuality.cpp" "mesh")
target_link_libraries(meshQualityCheck PRIVATE argparse::argparse nlohmann_json::nlohmann_json)
