set(AVAILABLE_LIBS animationIO)

if (NOT TARGET animationIO)
  message(STATUS "convertAnimation tool is not included. Missing: animationIO")
  return()
endif()

add_libpgo_tools(convertAnimation convertAnimation.cpp "${AVAILABLE_LIBS}")
target_link_libraries(convertAnimation PRIVATE argparse::argparse)