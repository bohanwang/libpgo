if(TARGET nlohmann_json::nlohmann_json)
  return()
endif()

message(STATUS "Loading nlohmann_json...")

include(FetchContent)
FetchContent_Declare(
  nlohmann_json
  URL https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.tar.gz
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES nlohmann_json
)
FetchContent_MakeAvailable(nlohmann_json)

message(STATUS "Done.")