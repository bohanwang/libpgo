function(pgo_should_strip_release_binary output_variable)
  set(should_strip OFF)

  if(NOT MSVC
     AND CMAKE_BUILD_TYPE STREQUAL "Release"
     AND NOT PGO_ENABLE_RELEASE_DEBUG_INFO)
    set(should_strip ON)
  endif()

  set(${output_variable} ${should_strip} PARENT_SCOPE)
endfunction()
