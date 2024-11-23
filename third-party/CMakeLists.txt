add_subdirectory(CCD)
add_subdirectory(ASA)
add_subdirectory(stb)
add_subdirectory(tiny_obj_loader)
add_subdirectory(tetgen)

if(MSVC)
  add_subdirectory(gmp-msvc)
  add_subdirectory(mpfr-msvc)
endif()

if(PGO_UI)
  add_subdirectory(imgui)
endif()

