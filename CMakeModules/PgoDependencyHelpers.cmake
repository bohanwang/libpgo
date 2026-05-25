include_guard(GLOBAL)

include(FetchContent)

function(pgo_dep_option name type value doc)
  set(${name} ${value} CACHE ${type} "${doc}" FORCE)
endfunction()

macro(pgo_fetch_make_available name)
  FetchContent_MakeAvailable(${name})
endmacro()

macro(pgo_fetch_populate_compat name reason)
  message(STATUS "Using FetchContent_Populate(${name}): ${reason}")
  FetchContent_GetProperties(${name})

  if(NOT ${name}_POPULATED)
    cmake_policy(PUSH)
    if(POLICY CMP0169)
      cmake_policy(SET CMP0169 OLD)
    endif()
    FetchContent_Populate(${name})
    cmake_policy(POP)
  endif()
endmacro()

function(pgo_apply_patch source_dir patch_file reason)
  if(NOT EXISTS "${patch_file}")
    message(FATAL_ERROR "Patch file does not exist: ${patch_file}")
  endif()

  find_package(Git QUIET)
  if(NOT GIT_FOUND)
    message(FATAL_ERROR "Git is required to apply patch: ${patch_file}")
  endif()

  get_filename_component(_pgo_patch_name "${patch_file}" NAME)
  get_filename_component(_pgo_patch_ceiling "${source_dir}" DIRECTORY)
  set(_pgo_patch_git
    "${CMAKE_COMMAND}" -E env "GIT_CEILING_DIRECTORIES=${_pgo_patch_ceiling}"
    "${GIT_EXECUTABLE}" -C "${source_dir}")

  execute_process(
    COMMAND ${_pgo_patch_git} apply --check "${patch_file}"
    RESULT_VARIABLE _pgo_patch_check_result
    OUTPUT_VARIABLE _pgo_patch_check_output
    ERROR_VARIABLE _pgo_patch_check_error
  )

  if(_pgo_patch_check_result EQUAL 0)
    message(STATUS "Applying patch ${_pgo_patch_name}: ${reason}")
    execute_process(
      COMMAND ${_pgo_patch_git} apply "${patch_file}"
      RESULT_VARIABLE _pgo_patch_apply_result
      OUTPUT_VARIABLE _pgo_patch_apply_output
      ERROR_VARIABLE _pgo_patch_apply_error
    )
    if(NOT _pgo_patch_apply_result EQUAL 0)
      message(FATAL_ERROR "Failed to apply ${patch_file}:\n${_pgo_patch_apply_error}")
    endif()
    return()
  endif()

  execute_process(
    COMMAND ${_pgo_patch_git} apply --reverse --check "${patch_file}"
    RESULT_VARIABLE _pgo_patch_reverse_check_result
    OUTPUT_VARIABLE _pgo_patch_reverse_check_output
    ERROR_VARIABLE _pgo_patch_reverse_check_error
  )

  if(_pgo_patch_reverse_check_result EQUAL 0)
    message(STATUS "Patch already applied ${_pgo_patch_name}: ${reason}")
  else()
    message(FATAL_ERROR
      "Patch ${patch_file} does not apply cleanly to ${source_dir}.\n"
      "Forward check:\n${_pgo_patch_check_error}\n"
      "Reverse check:\n${_pgo_patch_reverse_check_error}")
  endif()
endfunction()
