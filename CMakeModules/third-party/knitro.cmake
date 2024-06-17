if(TARGET Knitro::Knitro)
  return()
endif()

include(FindPackageHandleStandardArgs)

if(WIN32)
  if(KNITRO_C_INCLUDE_DIR AND KNITRO_LIBRARY_LIB AND KNITRO_LIBRARY_DLL)
  else()    
    find_path(KNITRO_C_INCLUDE_DIR
      NAMES "knitro.h"
      PATHS
      ${KNITRO_LIBRARY_HINT}
      ENV KNITRODIR
      PATH_SUFFIXES include
    )

    message(STATUS "Found Knitro include: ${KNITRO_C_INCLUDE_DIR}")

    find_library(KNITRO_LIBRARY_LIB
      NAMES "knitro1400.lib" "knitro1320.lib" "knitro1240.lib" "knitro1230.lib" "knitro1220.lib" "knitro1211.lib" "knitro1200.lib" "knitro1110.lib" 
      PATHS
      ENV KNITRODIR
      ${KNITRO_LIBRARY_HINT}
      PATH_SUFFIXES lib
    )

    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib" ".dll")

    message(STATUS "Found Knitro lib: ${KNITRO_LIBRARY_LIB}")

    find_library(KNITRO_LIBRARY_DLL
      NAMES "knitro1400.dll" "knitro1320.dll" "knitro1240.dll" "knitro1230.dll" "knitro1220.dll" "knitro1211.dll" "knitro1200.dll" "knitro1110.dll"  
      PATHS
      ENV KNITRODIR
      ${KNITRO_LIBRARY_HINT}
      PATH_SUFFIXES lib
    )

    SET(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")

    message(STATUS "Found Knitro dll: ${KNITRO_LIBRARY_DLL}")
  endif()

  if (KNITRO_C_INCLUDE_DIR AND KNITRO_LIBRARY_LIB AND KNITRO_LIBRARY_DLL)
    find_package_handle_standard_args(Knitro DEFAULT_MSG KNITRO_C_INCLUDE_DIR KNITRO_LIBRARY_LIB KNITRO_LIBRARY_DLL)
    mark_as_advanced(KNITRO_C_INCLUDE_DIR)
    mark_as_advanced(KNITRO_LIBRARY_LIB)
    mark_as_advanced(KNITRO_LIBRARY_DLL)

    add_library(KNITRO_LIB SHARED IMPORTED GLOBAL)
    set_target_properties(KNITRO_LIB PROPERTIES IMPORTED_LOCATION ${KNITRO_LIBRARY_DLL})
    set_target_properties(KNITRO_LIB PROPERTIES IMPORTED_IMPLIB ${KNITRO_LIBRARY_LIB})
    target_include_directories(KNITRO_LIB INTERFACE ${KNITRO_C_INCLUDE_DIR})
    # target_include_directories(KNITRO_LIB INTERFACE ${KNITRO_CPP_INCLUDE_DIR})
    add_library(Knitro::Knitro ALIAS KNITRO_LIB)
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  if (KNITRO_C_INCLUDE_DIR AND KNITRO_LIBRARY_LIB)
  else()    
    find_path(KNITRO_C_INCLUDE_DIR
      NAMES "knitro.h"
      PATHS
      ENV KNITRODIR      
      ${KNITRO_LIBRARY_HINT}
      PATH_SUFFIXES include
    )

    message(STATUS "Found Knitro HEADERS: ${KNITRO_C_INCLUDE_DIR}")

    find_library(KNITRO_LIBRARY_LIB
      NAMES knitro
      PATHS
      ENV KNITRODIR      
      ${KNITRO_LIBRARY_HINT}
      PATH_SUFFIXES lib
    )

    message(STATUS "Found Knitro lib: ${KNITRO_LIBRARY_LIB}")
  endif()

  if (KNITRO_C_INCLUDE_DIR AND KNITRO_LIBRARY_LIB)
    find_package_handle_standard_args(Knitro DEFAULT_MSG KNITRO_C_INCLUDE_DIR KNITRO_LIBRARY_LIB)
    mark_as_advanced(KNITRO_C_INCLUDE_DIR)
    mark_as_advanced(KNITRO_LIBRARY_LIB)

    add_library(KNITRO_LIB SHARED IMPORTED GLOBAL)
    set_target_properties(KNITRO_LIB PROPERTIES IMPORTED_LOCATION ${KNITRO_LIBRARY_LIB})
    target_include_directories(KNITRO_LIB INTERFACE ${KNITRO_C_INCLUDE_DIR})
    # target_include_directories(KNITRO_LIB INTERFACE ${KNITRO_CPP_INCLUDE_DIR})
    add_library(Knitro::Knitro ALIAS KNITRO_LIB)
  endif()
elseif(APPLE)
endif()