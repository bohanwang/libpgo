if(TARGET igl::core)
    return()
endif()

message(STATUS "Loading libigl...")

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG main
)

pgo_fetch_make_available(libigl)

message(STATUS "Done.")
