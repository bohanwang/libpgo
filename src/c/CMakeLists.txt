set(LIBPGO_C_HEADERS
  pgo_c.h
)

set(LIBPGO_C_SOURCES
  pgo_c.cpp
)

add_library(pgo_c SHARED ${LIBPGO_C_SOURCES} ${LIBPGO_C_HEADERS})

set_target_properties(pgo_c PROPERTIES
  VERSION 0.0.2
  OUTPUT_NAME pgo_c
  SOVERSION 0
  PUBLIC_HEADER "pgo_c.h"
)

target_include_directories(pgo_c PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(pgo_c PUBLIC
  $<BUILD_INTERFACE:LIBPGO_C_BUILD_DLL>
  $<INSTALL_INTERFACE:LIBPGO_C_LOAD_DLL>
)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    target_link_libraries(pgo_c PRIVATE ${lib})
  endif()
endforeach()

add_flag_poxis(pgo_c PRIVATE -fvisibility=hidden)
add_link_flag(pgo_c GNU All -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/pgo_c.version)
add_link_flag(pgo_c AppleClang All -Wl,-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/pgo_c.symobl.txt)

# pog_c install
include(CMakePackageConfigHelpers)

install(TARGETS pgo_c
  EXPORT pgo_cTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# create (temporary) export target file during build
export(EXPORT pgo_cTargets
  NAMESPACE pgo::
  FILE ${CMAKE_CURRENT_BINARY_DIR}/pgo_cTargets.cmake)

# install export target, config and version files for find_package
install(EXPORT pgo_cTargets
  NAMESPACE pgo::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pgo/")

if(PGO_ENABLE_PYTHON)
  add_library(pgo_c_static STATIC ${LIBPGO_C_SOURCES} ${LIBPGO_C_HEADERS})
  target_include_directories(pgo_c_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

  foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
    if(TARGET ${lib})
      target_link_libraries(pgo_c_static PUBLIC ${lib})
    endif()
  endforeach()
endif()