set(ASA_HEADERS 
  asa047.hpp
)

set(ASA_SRCS
  asa047.cpp
)

add_library(ASA_lib STATIC ${ASA_HEADERS} ${ASA_SRCS})
target_include_directories(ASA_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ASA_lib PUBLIC compilation_flag)

add_library(ASA::ASA ALIAS ASA_lib)
set_property(TARGET ASA_lib PROPERTY FOLDER "third-party")