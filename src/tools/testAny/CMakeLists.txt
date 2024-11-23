set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

list(APPEND AVAILABLE_LIBS pgo_c_static)
# add_libpgo_tools(testAny testAny.cpp "${AVAILABLE_LIBS}")

add_executable(testAny testAny.cpp)
target_include_directories(testAny PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../c/)

foreach(d ${AVAILABLE_LIBS})
  message(STATUS "testAny: ${d}")
  target_link_libraries(testAny PRIVATE ${d})
endforeach()

fix_mkl_tbb_debug(testAny)
fix_glut32(testAny)

target_link_libraries(testAny PRIVATE compilation_flag)

if(PGO_RELEASE_MODE_DEBUG)
  target_link_libraries(testAny PRIVATE compilation_flag_for_debug)
endif()