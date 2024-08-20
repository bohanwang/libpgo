FILE(READ "${CMAKE_CURRENT_SOURCE_DIR}/imgui.h" IMGUI_HEADER)
STRING(REGEX MATCH "#define[ \t]+IMGUI_VERSION[ \t]+\"([0-9]+\\.[0-9]+)\"" IMGUI_VER_MATCH "${IMGUI_HEADER}")
set(IMGUI_VER ${CMAKE_MATCH_1})

set(IMGUI_HEADERS
  imconfig.h
  imgui.h
  imgui_internal.h
  imstb_rectpack.h
  imstb_textedit.h
  imstb_truetype.h
)

set(IMGUI_SRC
  imgui.cpp
  imgui_demo.cpp
  imgui_draw.cpp
  imgui_tables.cpp
  imgui_widgets.cpp
)

add_library(imgui_lib STATIC ${IMGUI_HEADERS} ${IMGUI_SRC})
target_include_directories(imgui_lib PUBLIC ./)
target_link_libraries(imgui_lib PUBLIC compilation_flag)

add_library(imgui::imgui ALIAS imgui_lib)
set_property(TARGET imgui_lib PROPERTY FOLDER "third-party/imgui")

if(TARGET glfw)
  add_library(imgui_glfw3 STATIC ./backends/imgui_impl_glfw.cpp ./backends/imgui_impl_glfw.h)

  target_include_directories(imgui_glfw3 PUBLIC ./backends)
  target_link_libraries(imgui_glfw3 PUBLIC imgui::imgui glfw)

  add_library(imgui::glfw3 ALIAS imgui_glfw3)
  set_property(TARGET imgui_glfw3 PROPERTY FOLDER "third-party/imgui")
endif()

if(TARGET Vulkan::Vulkan)
  add_library(imgui_vulkan STATIC ./backends/imgui_impl_vulkan.cpp ./backends/imgui_impl_vulkan.h)

  target_include_directories(imgui_vulkan PUBLIC ./backends)
  target_link_libraries(imgui_vulkan PUBLIC imgui::imgui Vulkan::Vulkan)

  add_library(imgui::vulkan ALIAS imgui_vulkan)
  set_property(TARGET imgui_vulkan PROPERTY FOLDER "third-party/imgui")
endif()