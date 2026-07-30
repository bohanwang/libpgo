include(FetchContent)
include(GoogleTest)

set(INSTALL_GTEST OFF CACHE BOOL "Disable GoogleTest installation" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock build" FORCE)
set(gtest_force_shared_crt ON CACHE BOOL "Use shared CRT for GoogleTest on Windows" FORCE)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
  URL_HASH SHA256=f179ec217f9b3b3f3c6e8b02d3e7eda997b49e4ce26d6b235c9053bec9c0bf9f
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(googletest)

if(TARGET pgo_windows_runtime)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}")
  set(CMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE PRE_TEST)
  add_dependencies(gtest_main pgo_windows_runtime)
endif()

set(LIBPGO_TEST_FIXTURE_DIR "${CMAKE_CURRENT_BINARY_DIR}/fixtures")
set(LIBPGO_TEST_TET_BOX_DIR "${LIBPGO_TEST_FIXTURE_DIR}/tet-box")
set(LIBPGO_TEST_CUBIC_BOX_DIR "${LIBPGO_TEST_FIXTURE_DIR}/cubic-box")
set(LIBPGO_TEST_TET_BOX_VEG "${LIBPGO_TEST_TET_BOX_DIR}/box.veg")
set(LIBPGO_TEST_TET_BOX_OBJ "${LIBPGO_TEST_TET_BOX_DIR}/box.obj")
set(LIBPGO_TEST_CUBIC_BOX_VEG "${LIBPGO_TEST_CUBIC_BOX_DIR}/box.veg")
set(LIBPGO_TEST_CUBIC_BOX_OBJ "${LIBPGO_TEST_CUBIC_BOX_DIR}/box.obj")
set(LIBPGO_TEST_TET_BOX_JSON "${LIBPGO_TEST_TET_BOX_DIR}/box.json")
set(LIBPGO_TEST_CUBIC_BOX_JSON "${LIBPGO_TEST_CUBIC_BOX_DIR}/box.json")

file(MAKE_DIRECTORY "${LIBPGO_TEST_TET_BOX_DIR}" "${LIBPGO_TEST_CUBIC_BOX_DIR}")
configure_file(
  "${CMAKE_SOURCE_DIR}/examples/assets/volume/box/box-tet.veg"
  "${LIBPGO_TEST_TET_BOX_VEG}"
  COPYONLY)
configure_file(
  "${CMAKE_SOURCE_DIR}/examples/assets/volume/box/box.obj"
  "${LIBPGO_TEST_TET_BOX_OBJ}"
  COPYONLY)
configure_file(
  "${CMAKE_SOURCE_DIR}/examples/assets/volume/box/box.obj"
  "${LIBPGO_TEST_CUBIC_BOX_OBJ}"
  COPYONLY)
configure_file(
  "${CMAKE_SOURCE_DIR}/examples/assets/volume/box/box-fixed.txt"
  "${LIBPGO_TEST_TET_BOX_DIR}/box-fixed.txt"
  COPYONLY)
configure_file(
  "${CMAKE_SOURCE_DIR}/examples/assets/volume/box/box-fixed.txt"
  "${LIBPGO_TEST_CUBIC_BOX_DIR}/box-fixed.txt"
  COPYONLY)

function(libpgo_materialize_box_config source destination)
  file(READ "${source}" config_contents)
  string(REGEX REPLACE
    "\"(tet|cubic)-mesh\"[ \t]*:[ \t]*\"[^\"]+\""
    "\"\\1-mesh\": \"box.veg\""
    config_contents "${config_contents}")
  string(REGEX REPLACE
    "\"surface-mesh\"[ \t]*:[ \t]*\"[^\"]+\""
    "\"surface-mesh\": \"box.obj\""
    config_contents "${config_contents}")
  string(REGEX REPLACE
    "\"filename\"[ \t]*:[ \t]*\"[^\"]*box-fixed.txt\""
    "\"filename\": \"box-fixed.txt\""
    config_contents "${config_contents}")
  file(WRITE "${destination}" "${config_contents}")
endfunction()

libpgo_materialize_box_config(
  "${CMAKE_SOURCE_DIR}/examples/configs/volume/box/box-tet-sampled.json"
  "${LIBPGO_TEST_TET_BOX_JSON}")
libpgo_materialize_box_config(
  "${CMAKE_SOURCE_DIR}/examples/configs/volume/box/box-cubic-sampled.json"
  "${LIBPGO_TEST_CUBIC_BOX_JSON}")
libpgo_materialize_box_config(
  "${CMAKE_SOURCE_DIR}/examples/configs/volume/box/box-hang-tet-ipc.json"
  "${LIBPGO_TEST_TET_BOX_DIR}/box-ipc.json")
libpgo_materialize_box_config(
  "${CMAKE_SOURCE_DIR}/examples/configs/volume/box/box-hang-cubic-ipc.json"
  "${LIBPGO_TEST_CUBIC_BOX_DIR}/box-ipc.json")

add_custom_command(
  OUTPUT "${LIBPGO_TEST_CUBIC_BOX_VEG}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${LIBPGO_TEST_CUBIC_BOX_DIR}"
  COMMAND $<TARGET_FILE:cubicMesher>
    --input-mesh "${LIBPGO_TEST_CUBIC_BOX_OBJ}"
    --resolution 4
    --output-mesh "${LIBPGO_TEST_CUBIC_BOX_VEG}"
    --E 10000000
    --nu 0.45
    --density 1000
  DEPENDS
    cubicMesher
    "${LIBPGO_TEST_CUBIC_BOX_OBJ}"
  VERBATIM)
add_custom_target(libpgo_test_cubic_assets
  DEPENDS "${LIBPGO_TEST_CUBIC_BOX_VEG}")

add_executable(pgo_c_gtest pgo_c_gtest.cpp)
target_include_directories(pgo_c_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/c)
target_link_libraries(pgo_c_gtest PRIVATE GTest::gtest_main pgo_c)
target_compile_definitions(pgo_c_gtest PRIVATE
  LIBPGO_TEST_TORUS_VEG="${CMAKE_SOURCE_DIR}/examples/assets/common/torus-tet.veg")
set_property(TARGET pgo_c_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(pgo_c_gtest)

add_executable(mesh_gtest mesh_gtest.cpp)
target_include_directories(mesh_gtest PRIVATE ${CMAKE_SOURCE_DIR}/src/core/mesh)
target_link_libraries(mesh_gtest PRIVATE GTest::gtest_main mesh)
set_property(TARGET mesh_gtest PROPERTY FOLDER "tests/gtest")
gtest_discover_tests(mesh_gtest)

add_subdirectory(src)
