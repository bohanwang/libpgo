set(PGO_PYTHON_SOURCES pypgo CACHE INTERNAL "")
set(PGO_PYLIST_FILE "${CMAKE_CURRENT_BINARY_DIR}/pypgo_list.cpp" CACHE INTERNAL "")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/pypgo_list.cpp" "")

add_library(pypgo MODULE pypgo.cpp pypgo.h pypgo_main.cpp "${CMAKE_CURRENT_BINARY_DIR}/pypgo_list.cpp")
target_include_directories(pypgo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pypgo PRIVATE pybind11::module pybind11::lto pybind11::windows_extras)
target_link_libraries(pypgo PRIVATE pgo_c_static)
target_compile_definitions(pypgo PRIVATE VERSION_INFO=${PYPGO_VERSION_INFO})

pybind11_extension(pypgo)

if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS
    pybind11_strip(pypgo)
endif()

set_target_properties(pypgo PROPERTIES CXX_VISIBILITY_PRESET "hidden"
    CUDA_VISIBILITY_PRESET "hidden")

if(TARGET pgo_windows_runtime)
    add_dependencies(pypgo pgo_windows_runtime)
    add_custom_command(
        TARGET pypgo
        POST_BUILD
        COMMAND
            ${CMAKE_COMMAND} -E copy_directory_if_different
            "${PGO_WINDOWS_RUNTIME_OUTPUT_DIRECTORY}"
            "$<TARGET_FILE_DIR:pypgo>"
        COMMENT "Staging Windows runtime DLLs beside $<TARGET_FILE_NAME:pypgo>"
        VERBATIM)
endif()

if(PGO_PORTABLE_BUILD)
    if(APPLE)
        set(PGO_PYTHON_RUNTIME_PATH "@loader_path/../..")
    elseif(UNIX)
        set(PGO_PYTHON_RUNTIME_PATH "$ORIGIN/../..")
    endif()

    if(DEFINED PGO_PYTHON_RUNTIME_PATH)
        set_target_properties(pypgo PROPERTIES
            BUILD_RPATH "${PGO_PYTHON_RUNTIME_PATH}"
            BUILD_WITH_INSTALL_RPATH TRUE
            INSTALL_RPATH "${PGO_PYTHON_RUNTIME_PATH}")
    endif()
endif()
