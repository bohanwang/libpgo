if(TARGET Pardiso::Pardiso)
  return()
endif()

include(FindPackageHandleStandardArgs)

if(WIN32)  
  message(STATUS "PARDISO support is not available on Windows")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  if (PARDISO_C_INCLUDE_DIR AND PARDISO_LIBRARY_LIB)
  else()    
    find_path(PARDISO_C_INCLUDE_DIR
      NAMES "pardiso.h"
      PATHS
      ENV PARDISOROOT      
      ${PARDISO_LIBRARY_HINT}
      PATH_SUFFIXES include
    )

    message(STATUS "Found PARDISO HEADERS: ${PARDISO_C_INCLUDE_DIR}")

    find_library(PARDISO_LIBRARY_LIB
      NAMES pardiso
      PATHS
      ENV PARDISOROOT      
      ${PARDISO_LIBRARY_HINT}
      PATH_SUFFIXES lib
    )

    message(STATUS "Found PARDISO lib: ${PARDISO_LIBRARY_LIB}")
  endif()

  if (PARDISO_C_INCLUDE_DIR AND PARDISO_LIBRARY_LIB)
    find_package_handle_standard_args(pardiso DEFAULT_MSG PARDISO_C_INCLUDE_DIR PARDISO_LIBRARY_LIB)
    mark_as_advanced(PARDISO_C_INCLUDE_DIR)
    mark_as_advanced(PARDISO_LIBRARY_LIB)

    add_library(PARDISO_LIB SHARED IMPORTED GLOBAL)
    set_target_properties(PARDISO_LIB PROPERTIES IMPORTED_LOCATION ${PARDISO_LIBRARY_LIB})
    target_include_directories(PARDISO_LIB INTERFACE ${PARDISO_C_INCLUDE_DIR})
    # target_include_directories(PARDISO_LIB INTERFACE ${PARDISO_CPP_INCLUDE_DIR})

    target_link_libraries(PARDISO_LIB INTERFACE MKL::MKL Ceres::ceres)
    add_library(Pardiso::Pardiso ALIAS PARDISO_LIB)
  endif()
elseif(APPLE)
  message(STATUS "PARDISO support is not available on macOS")
endif()