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

FetchContent_MakeAvailable(libigl)

message(STATUS "Done.")