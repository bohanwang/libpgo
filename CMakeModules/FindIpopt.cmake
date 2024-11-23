include(FindPackageHandleStandardArgs)

# IpTNLP.hpp

if(WIN32)
  if (IPOPT_LIBRARY_LIB_DEBUG AND
    IPOPT_LIBRARY_LIB_RELEASE AND
    IPOPT_INCLUDE)
  else()
    find_library(IPOPT_LIBRARY_LIB_DEBUG
      NAMES "libipopt" "ipopt.lib"
      PATHS
      ${EXT_LIB_ROOT}
      ${EXT_LIB_ROOT}/lib
      ${IPOPT_LIBRARY_HINT}
      ${IPOPT_LIBRARY_HINT}/lib
      PATH_SUFFIXES Debug
    )

    find_library(IPOPT_LIBRARY_LIB_RELEASE
      NAMES "libipopt" "ipopt.lib"
      PATHS
      ${EXT_LIB_ROOT}
      ${EXT_LIB_ROOT}/lib
      ${IPOPT_LIBRARY_HINT}
      ${IPOPT_LIBRARY_HINT}/lib
      PATH_SUFFIXES Release
    )

    find_path(IPOPT_INCLUDE
      NAMES "coin-or/IpTNLP.hpp" 
      PATHS
      ${EXT_LIB_ROOT}
      ${EXT_LIB_ROOT}/include
      ${IPOPT_LIBRARY_HINT}
      ${IPOPT_LIBRARY_HINT}/include
    )

  endif()

  if (IPOPT_LIBRARY_LIB_DEBUG AND
    IPOPT_LIBRARY_LIB_RELEASE AND
    IPOPT_INCLUDE)

    find_package_handle_standard_args(Ipopt DEFAULT_MSG IPOPT_LIBRARY_LIB_DEBUG IPOPT_LIBRARY_LIB_RELEASE IPOPT_INCLUDE)
    mark_as_advanced(IPOPT_LIBRARY_LIB_DEBUG)
    mark_as_advanced(IPOPT_LIBRARY_LIB_RELEASE)
    mark_as_advanced(IPOPT_INCLUDE)

    add_library(Ipopt::Core STATIC IMPORTED GLOBAL)
    target_include_directories(Ipopt::Core INTERFACE ${IPOPT_INCLUDE})
    set_target_properties(Ipopt::Core PROPERTIES IMPORTED_LOCATION_DEBUG ${IPOPT_LIBRARY_LIB_DEBUG})
    set_target_properties(Ipopt::Core PROPERTIES IMPORTED_LOCATION_RELEASE ${IPOPT_LIBRARY_LIB_RELEASE})

    target_link_libraries(Ipopt::Core INTERFACE MKL::MKL)
  endif()
else()
  if (IPOPT_LIBRARY_LIB AND
    IPOPT_INCLUDE)
  else()
    find_library(IPOPT_LIBRARY_LIB
      NAMES "ipopt"
      PATHS
      ${EXT_LIB_ROOT}
      ${EXT_LIB_ROOT}/lib
      ${IPOPT_LIBRARY_HINT}
      ${IPOPT_LIBRARY_HINT}/lib
    )

    find_path(IPOPT_INCLUDE
      NAMES "coin-or/IpTNLP.hpp"
      PATHS
      ${EXT_LIB_ROOT}
      ${EXT_LIB_ROOT}/include
      ${IPOPT_LIBRARY_HINT}
      ${IPOPT_LIBRARY_HINT}/include
    )
  endif()

  if (IPOPT_LIBRARY_LIB AND
    IPOPT_INCLUDE)

    find_package_handle_standard_args(Ipopt DEFAULT_MSG IPOPT_LIBRARY_LIB IPOPT_INCLUDE)
    mark_as_advanced(IPOPT_LIBRARY_LIB)
    mark_as_advanced(IPOPT_INCLUDE)

    add_library(Ipopt::Core STATIC IMPORTED GLOBAL)
    target_include_directories(Ipopt::Core INTERFACE ${IPOPT_INCLUDE})
    set_target_properties(Ipopt::Core PROPERTIES IMPORTED_LOCATION ${IPOPT_LIBRARY_LIB})
    target_link_libraries(Ipopt::Core INTERFACE MKL::MKL)
  endif()
endif()
