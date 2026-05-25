if(TARGET glfw)
  return()
endif()

message(STATUS "Loading glfw...")

pgo_dep_option(BUILD_SHARED_LIBS BOOL OFF "Build shared libraries")

pgo_dep_option(GLFW_BUILD_EXAMPLES BOOL OFF "Build the GLFW example programs")
pgo_dep_option(GLFW_BUILD_TESTS BOOL OFF "Build the GLFW test programs")
pgo_dep_option(GLFW_BUILD_DOCS BOOL OFF "Build the GLFW documentation")
pgo_dep_option(GLFW_INSTALL BOOL OFF "Generate installation target")

include(FetchContent)
FetchContent_Declare(
  glfw
  URL https://github.com/glfw/glfw/releases/download/3.4/glfw-3.4.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(glfw)

message(STATUS "Done.")
