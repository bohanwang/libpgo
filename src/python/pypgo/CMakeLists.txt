set(PGO_PYTHON_SOURCES pypgo CACHE INTERNAL "")
set(PGO_PYLIST_FILE "${CMAKE_CURRENT_BINARY_DIR}/pypgo_list.cpp" CACHE INTERNAL "")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pypgo_list.cpp" "")

add_library(pypgo MODULE pypgo.cpp pypgo.h pypgo_main.cpp "${CMAKE_CURRENT_BINARY_DIR}/pypgo_list.cpp")
target_include_directories(pypgo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pypgo PRIVATE pybind11::module pybind11::lto pybind11::windows_extras)
target_link_libraries(pypgo PRIVATE pgo_c_static)
target_compile_definitions(pypgo PRIVATE VERSION_INFO=${PYPGO_VERSION_INFO})

pybind11_extension(pypgo)

pgo_should_strip_release_binary(PGO_STRIP_PYTHON_EXTENSION)
if(PGO_STRIP_PYTHON_EXTENSION)
    # Keep Release symbols for source diagnostics, but strip release wheels.
    pybind11_strip(pypgo)
endif()

set_target_properties(pypgo PROPERTIES OUTPUT_NAME "_pypgo"
    CXX_VISIBILITY_PRESET "hidden"
    CUDA_VISIBILITY_PRESET "hidden")

add_custom_command(
    TARGET pypgo
    POST_BUILD
    COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/pgo_run_sim.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/pgo_dump_abc.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/pgo_run_cases.py"
        "$<TARGET_FILE_DIR:pypgo>"
    COMMENT "Staging the pypgo Python package beside $<TARGET_FILE_NAME:pypgo>"
    VERBATIM)

pgo_attach_runtime_dependencies(pypgo COPY_RUNTIME WHEEL_EXTENSION)
