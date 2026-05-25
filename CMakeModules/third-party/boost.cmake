if(TARGET Boost::boost)
  return()
endif()

message(STATUS "Loading Boost...")
include(FetchContent)

pgo_dep_option(BOOST_INCLUDE_LIBRARIES STRING
  "any;foreach;format;graph;heap;logic;math;multiprecision;property_map;system;thread;variant"
  "Boost libraries used by CGAL/libpgo")

# set(FETCHCONTENT_QUIET OFF)
FetchContent_Declare(
  boost
  URL https://github.com/boostorg/boost/releases/download/boost-1.85.0/boost-1.85.0-cmake.tar.xz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(boost)

if(DEFINED boost_SOURCE_DIR)
  pgo_dep_option(Boost_NO_SYSTEM_PATHS BOOL ON "Restrict Boost lookup to fetched Boost")
  pgo_dep_option(Boost_NO_BOOST_CMAKE BOOL ON "Avoid external boost-cmake package lookup")
  pgo_dep_option(Boost_INCLUDE_DIR PATH "${boost_SOURCE_DIR}/libs/config/include" "Fetched Boost include directory for FindBoost compatibility")
  pgo_dep_option(Boost_INCLUDE_DIRS STRING "${boost_SOURCE_DIR}/libs/config/include" "Fetched Boost include directories for FindBoost compatibility")
endif()

message(STATUS "Done.")
