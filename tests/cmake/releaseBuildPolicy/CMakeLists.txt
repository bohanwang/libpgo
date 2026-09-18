cmake_minimum_required(VERSION 3.29)
project(releaseBuildPolicySmoke LANGUAGES NONE)

if(NOT DEFINED LIBPGO_SOURCE_DIR)
  message(FATAL_ERROR "LIBPGO_SOURCE_DIR is required")
endif()

if(NOT DEFINED EXPECTED_STRIP)
  message(FATAL_ERROR "EXPECTED_STRIP is required")
endif()

include("${LIBPGO_SOURCE_DIR}/CMakeModules/releaseBuildPolicy.cmake")
pgo_should_strip_release_binary(actual_strip)

if(NOT "${actual_strip}" STREQUAL "${EXPECTED_STRIP}")
  message(FATAL_ERROR
    "Expected strip policy ${EXPECTED_STRIP}, got ${actual_strip}")
endif()
