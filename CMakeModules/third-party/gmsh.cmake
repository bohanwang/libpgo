if(TARGET Gmsh::Gmsh)
    return()
endif()

include(FindPackageHandleStandardArgs)

set(_pgo_gmsh_prefix_hints)

if(PGO_CHECK_CONDA AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
    if(WIN32)
        list(APPEND _pgo_gmsh_prefix_hints "$ENV{CONDA_PREFIX}/Library")
    else()
        list(APPEND _pgo_gmsh_prefix_hints "$ENV{CONDA_PREFIX}")
    endif()
endif()

list(APPEND _pgo_gmsh_prefix_hints
    ${GMSH_ROOT}
    $ENV{GMSH_ROOT}
    ${GMSH_LIBRARY_HINT}
)
list(REMOVE_DUPLICATES _pgo_gmsh_prefix_hints)

function(_pgo_create_gmsh_alias target_name)
    if(TARGET Gmsh::Gmsh)
        return()
    endif()

    add_library(Gmsh::Gmsh INTERFACE IMPORTED GLOBAL)
    target_link_libraries(Gmsh::Gmsh INTERFACE ${target_name})
endfunction()

function(_pgo_alias_existing_gmsh_target result_var)
    foreach(_pgo_gmsh_target IN ITEMS Gmsh::Gmsh gmsh::shared gmsh::lib gmsh)
        if(TARGET ${_pgo_gmsh_target})
            if(NOT _pgo_gmsh_target STREQUAL "Gmsh::Gmsh")
                _pgo_create_gmsh_alias(${_pgo_gmsh_target})
            endif()
            set(${result_var} TRUE PARENT_SCOPE)
            return()
        endif()
    endforeach()

    set(${result_var} FALSE PARENT_SCOPE)
endfunction()

find_package(Gmsh CONFIG QUIET
    PATHS ${_pgo_gmsh_prefix_hints}
)
_pgo_alias_existing_gmsh_target(_pgo_gmsh_found_config_target)
if(_pgo_gmsh_found_config_target)
    return()
endif()

find_package(gmsh CONFIG QUIET
    PATHS ${_pgo_gmsh_prefix_hints}
)
_pgo_alias_existing_gmsh_target(_pgo_gmsh_found_config_target)
if(_pgo_gmsh_found_config_target)
    return()
endif()

if(NOT GMSH_INCLUDE_DIR)
    find_path(GMSH_INCLUDE_DIR
        NAMES "gmsh.h"
        PATHS ${_pgo_gmsh_prefix_hints}
        PATH_SUFFIXES include
    )

    message(STATUS "Found GMSH HEADERS: ${GMSH_INCLUDE_DIR}")
endif()

if(WIN32)
    if(NOT GMSH_LIBRARY)
        find_library(GMSH_LIBRARY
            NAMES gmsh
            PATHS ${_pgo_gmsh_prefix_hints}
            PATH_SUFFIXES lib
        )

        message(STATUS "Found GMSH import lib: ${GMSH_LIBRARY}")
    endif()

    if(NOT GMSH_RUNTIME_LIBRARY)
        find_file(GMSH_RUNTIME_LIBRARY
            NAMES gmsh.dll
            PATHS ${_pgo_gmsh_prefix_hints}
            PATH_SUFFIXES bin
        )

        message(STATUS "Found GMSH runtime library: ${GMSH_RUNTIME_LIBRARY}")
    endif()

    if(GMSH_INCLUDE_DIR AND GMSH_LIBRARY AND GMSH_RUNTIME_LIBRARY)
        find_package_handle_standard_args(GMSH DEFAULT_MSG
            GMSH_INCLUDE_DIR
            GMSH_LIBRARY
            GMSH_RUNTIME_LIBRARY
        )
        mark_as_advanced(GMSH_INCLUDE_DIR)
        mark_as_advanced(GMSH_LIBRARY)
        mark_as_advanced(GMSH_RUNTIME_LIBRARY)

        add_library(GMSH_LIB SHARED IMPORTED GLOBAL)
        set_target_properties(GMSH_LIB PROPERTIES
            IMPORTED_IMPLIB ${GMSH_LIBRARY}
            IMPORTED_LOCATION ${GMSH_RUNTIME_LIBRARY}
        )
        target_include_directories(GMSH_LIB INTERFACE ${GMSH_INCLUDE_DIR})
        add_library(Gmsh::Gmsh ALIAS GMSH_LIB)
    endif()
else()
    if(NOT GMSH_LIBRARY)
        find_library(GMSH_LIBRARY
            NAMES gmsh
            PATHS ${_pgo_gmsh_prefix_hints}
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
endif()
