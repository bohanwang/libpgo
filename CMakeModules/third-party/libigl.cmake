if(TARGET igl::core)
    return()
endif()

message(STATUS "Loading libigl...")

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG 477e15a3d566a21f415aa5ee62992b12a836b01b
)

FetchContent_MakeAvailable(libigl)

message(STATUS "Done.")
