if(TARGET Boost::thread AND TARGET Boost::system AND TARGET Boost::iostreams)
  return()
endif()

message(STATUS "Loading Boost from conda...")

if(NOT PGO_CHECK_CONDA OR "$ENV{CONDA_PREFIX}" STREQUAL "")
  message(FATAL_ERROR "Boost is required from the active conda environment. Activate conda and install libboost-devel.")
endif()

if(WIN32)
  set(PGO_BOOST_PREFIX "$ENV{CONDA_PREFIX}/Library")
else()
  set(PGO_BOOST_PREFIX "$ENV{CONDA_PREFIX}")
endif()

list(PREPEND CMAKE_PREFIX_PATH "${PGO_BOOST_PREFIX}")

set(BOOST_ROOT "${PGO_BOOST_PREFIX}" CACHE PATH "Boost prefix" FORCE)
set(BOOST_INCLUDEDIR "${PGO_BOOST_PREFIX}/include" CACHE PATH "Boost include directory" FORCE)
set(BOOST_LIBRARYDIR "${PGO_BOOST_PREFIX}/lib" CACHE PATH "Boost library directory" FORCE)
set(Boost_NO_SYSTEM_PATHS ON CACHE BOOL "Restrict Boost lookup to conda" FORCE)
set(Boost_NO_BOOST_CMAKE OFF CACHE BOOL "Prefer conda Boost CMake config" FORCE)

find_package(Boost CONFIG REQUIRED COMPONENTS
  any
  foreach
  format
  graph
  heap
  iostreams
  logic
  math
  multiprecision
  property_map
  system
  thread
  variant
)

message(STATUS "Done.")
