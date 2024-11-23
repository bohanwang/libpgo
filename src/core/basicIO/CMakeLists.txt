if(NOT TARGET pgoLogging)
  message(STATUS "basicIO is not included. Missing: pgoLogging")
  return()
endif()

add_libpgo_lib(basicIO basicIO.cpp basicIO.h)
target_link_libraries(basicIO PUBLIC pgoLogging)