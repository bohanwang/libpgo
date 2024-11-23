project(CCD LANGUAGES CXX C)

# Safe CCD
add_library(CCD_SafeCCD STATIC ./SafeCCD/SAFE_CCD.h ./SafeCCD/MY_MATH.h ./SafeCCD/safeCCDWrapper.h ./SafeCCD/safeCCDWrapper.cpp)
target_include_directories(CCD_SafeCCD PUBLIC ./SafeCCD)
target_link_libraries(CCD_SafeCCD PUBLIC compilation_flag)
set_property(TARGET CCD_SafeCCD PROPERTY CXX_STANDARD 17)

add_library(CCD::SafeCCD ALIAS CCD_SafeCCD)
set_property(TARGET CCD_SafeCCD PROPERTY FOLDER "third-party/CCD")

# Eaxct CCD
add_library(CCD_ExactCCD STATIC ./ExactCCD/expansion.h ./ExactCCD/interval.h ./ExactCCD/vec.h ./ExactCCD/rootparitycollisiontest.h ./ExactCCD/interval_base.h 
  ./ExactCCD/expansion.cpp ./ExactCCD/interval.cpp ./ExactCCD/rootparitycollisiontest.cpp
  ./ExactCCD/exactCCDWrapper.h ./ExactCCD/exactCCDWrapper.cpp)

add_library(CCD::Exact ALIAS CCD_ExactCCD)

target_link_libraries(CCD_ExactCCD PUBLIC compilation_flag)
target_include_directories(CCD_ExactCCD PUBLIC ./ExactCCD)

set_property(TARGET CCD_ExactCCD PROPERTY FOLDER "third-party/CCD")
