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