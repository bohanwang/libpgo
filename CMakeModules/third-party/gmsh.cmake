if(TARGET Gmsh::Gmsh)
    return()
endif()

include(FindPackageHandleStandardArgs)

if(WIN32)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(GMSH_INCLUDE_DIR AND GMSH_LIBRARY)
    else()
        find_path(GMSH_INCLUDE_DIR
            NAMES "gmsh.h"
            PATH_SUFFIXES include
        )

        message(STATUS "Found GMSH HEADERS: ${GMSH_INCLUDE_DIR}")

        find_library(GMSH_LIBRARY
            NAMES gmsh
            PATHS
            PATH_SUFFIXES lib
        )

        message(STATUS "Found GMSH lib: ${GMSH_LIBRARY}")
    endif()

    if(GMSH_INCLUDE_DIR AND GMSH_LIBRARY)
        find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH_INCLUDE_DIR GMSH_LIBRARY)
        mark_as_advanced(GMSH_INCLUDE_DIR)
        mark_as_advanced(GMSH_LIBRARY)

        add_library(GMSH_LIB SHARED IMPORTED GLOBAL)
        set_target_properties(GMSH_LIB PROPERTIES IMPORTED_LOCATION ${GMSH_LIBRARY})
        target_include_directories(GMSH_LIB INTERFACE ${GMSH_INCLUDE_DIR})
        add_library(Gmsh::Gmsh ALIAS GMSH_LIB)
    endif()
elseif(APPLE)
endif()
