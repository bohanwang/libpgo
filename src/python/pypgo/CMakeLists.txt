

add_library(pypgo MODULE src/pypgo.cpp)
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