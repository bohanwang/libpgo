if(TARGET nanobind)
  return()
endif()

message(STATUS "Downloading nanobind...")

include(FetchContent)
FetchContent_Declare(
  nanobind
  GIT_REPOSITORY https://github.com/wjakob/nanobind.git
  GIT_TAG v2.12.0
  GIT_SUBMODULES ext/robin_map
  GIT_SHALLOW TRUE
)

pgo_fetch_make_available(nanobind)

message(STATUS "nanobind ready.")
