if(TARGET Boost::boost)
  return()
endif()

message(STATUS "Loading Boost...")
include(FetchContent)

# set(FETCHCONTENT_QUIET OFF)
FetchContent_Declare(
  boost
  URL https://github.com/boostorg/boost/releases/download/boost-1.85.0/boost-1.85.0-cmake.tar.xz
  URL_HASH SHA256=0a9cc56ceae46986f5f4d43fe0311d90cf6d2fa9028258a95cab49ffdacf92ad
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  OVERRIDE_FIND_PACKAGE
  # FIND_PACKAGE_ARGS NAMES Boost COMPONENTS any foreach format graph heap logic math multiprecision property_map system thread variant system program_options serialization
)

FetchContent_MakeAvailable(boost)

message(STATUS "Done.")
