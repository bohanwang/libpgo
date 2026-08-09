function(pgo_read_project_version pyproject_path output_variable)
  file(STRINGS "${pyproject_path}"
    _version_lines
    REGEX "^version[ \t]*=[ \t]*\"[0-9]+\\.[0-9]+\\.[0-9]+(\\.[0-9]+)?\"[ \t]*$")

  list(LENGTH _version_lines _version_count)
  if(NOT _version_count EQUAL 1)
    message(FATAL_ERROR
      "Expected exactly one CMake-compatible [project].version in ${pyproject_path}")
  endif()

  list(GET _version_lines 0 _version_line)
  string(REGEX REPLACE
    "^version[ \t]*=[ \t]*\"([0-9]+\\.[0-9]+\\.[0-9]+(\\.[0-9]+)?)\"[ \t]*$"
    "\\1"
    _version
    "${_version_line}")

  set("${output_variable}" "${_version}" PARENT_SCOPE)
endfunction()
