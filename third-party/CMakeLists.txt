if(NOT TARGET glad)
  add_subdirectory(glad)
endif()

# ##############################
# # libigl
# set(LIBIGL_USE_STATIC_LIBRARY OFF CACHE BOOL "Use libigl as static library" FORCE)
# set(LIBIGL_COPYLEFT_CGAL ON CACHE BOOL "Use CGAL" FORCE)
# set(LIBIGL_COPYLEFT_COMISO OFF CACHE BOOL "Use CoMiso" FORCE)
# set(LIBIGL_EMBREE OFF CACHE BOOL "Use Embree" FORCE)
# set(LIBIGL_RESTRICTED_MATLAB OFF CACHE BOOL "Use Matlab" FORCE)
# set(LIBIGL_RESTRICTED_MOSEK OFF CACHE BOOL "Use MOSEK" FORCE)
# set(LIBIGL_OPENGL OFF CACHE BOOL "Use OpenGL" FORCE)
# set(LIBIGL_GLFW OFF CACHE BOOL "Use GLFW" FORCE)
# set(LIBIGL_IMGUI OFF CACHE BOOL "Use ImGui" FORCE)
# set(LIBIGL_PNG OFF CACHE BOOL "Use PNG" FORCE)
# set(LIBIGL_COPYLEFT_TETGEN OFF CACHE BOOL "Use Tetgen" FORCE)
# set(LIBIGL_RESTRICTED_TRIANGLE OFF CACHE BOOL "Use Triangle" FORCE)
# set(LIBIGL_PREDICATES OFF CACHE BOOL "Use exact predicates" FORCE)
# set(LIBIGL_WITH_XML OFF CACHE BOOL "Use XML" FORCE)
# set(LIBIGL_INSTALL OFF CACHE BOOL "Export libigl CMake targets" FORCE)

# set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/libigl/cmake")

# add_subdirectory(tetgen)

# add_subdirectory(libigl)

add_subdirectory(CCD)

# add_subdirectory(libShapeOp)

# add_subdirectory(autodiff/autodiff)

# add_subdirectory(ASA)

add_subdirectory(tiny_obj_loader)

# add_subdirectory(stb)

# add_library(GTE INTERFACE)
# target_include_directories(GTE INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/GTE)

# set(TETWILD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/TetWild/extern")
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "TetWild/cmake")

# # Generate position independent code
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ################################################################################
# # 3rd party libraries
# ################################################################################
# if(NOT TARGET geogram)
# tetwild_download_geogram()
# include(geogram)
# endif()

# add_subdirectory(${TETWILD_EXTERNAL}/pymesh)

# ################################################################################
# # TetWild
# ################################################################################
# # Build static library for executable
# add_library(libTetWild
# TetWild/include/tetwild/Args.h
# TetWild/include/tetwild/Exception.h
# TetWild/include/tetwild/Logger.h
# TetWild/include/tetwild/tetwild.h
# TetWild/src/tetwild/BSPSubdivision.cpp
# TetWild/src/tetwild/BSPSubdivision.h
# TetWild/src/tetwild/CGALTypes.h
# TetWild/src/tetwild/Common.cpp
# TetWild/src/tetwild/Common.h
# TetWild/src/tetwild/DelaunayTetrahedralization.cpp
# TetWild/src/tetwild/DelaunayTetrahedralization.h
# TetWild/src/tetwild/DistanceQuery.cpp
# TetWild/src/tetwild/DistanceQuery.h
# TetWild/src/tetwild/EdgeCollapser.cpp
# TetWild/src/tetwild/EdgeCollapser.h
# TetWild/src/tetwild/EdgeRemover.cpp
# TetWild/src/tetwild/EdgeRemover.h
# TetWild/src/tetwild/EdgeSplitter.cpp
# TetWild/src/tetwild/EdgeSplitter.h
# TetWild/src/tetwild/ForwardDecls.h
# TetWild/src/tetwild/InoutFiltering.cpp
# TetWild/src/tetwild/InoutFiltering.h
# TetWild/src/tetwild/LocalOperations.cpp
# TetWild/src/tetwild/LocalOperations.h
# TetWild/src/tetwild/Logger.cpp
# TetWild/src/tetwild/MeshConformer.cpp
# TetWild/src/tetwild/MeshConformer.h
# TetWild/src/tetwild/MeshRefinement.cpp
# TetWild/src/tetwild/MeshRefinement.h
# TetWild/src/tetwild/Preprocess.cpp
# TetWild/src/tetwild/Preprocess.h
# TetWild/src/tetwild/SimpleTetrahedralization.cpp
# TetWild/src/tetwild/SimpleTetrahedralization.h
# TetWild/src/tetwild/State.cpp
# TetWild/src/tetwild/State.h
# TetWild/src/tetwild/TetmeshElements.cpp
# TetWild/src/tetwild/TetmeshElements.h
# TetWild/src/tetwild/tetwild.cpp
# TetWild/src/tetwild/VertexSmoother.cpp
# TetWild/src/tetwild/VertexSmoother.h
# TetWild/src/tetwild/geogram/mesh_AABB.cpp
# TetWild/src/tetwild/geogram/mesh_AABB.h
# )

# target_include_directories(libTetWild
# PRIVATE
# TetWild/src
# PUBLIC
# TetWild/include
# )

# target_link_libraries(libTetWild
# PUBLIC
# geogram
# igl::core
# pymesh::pymesh
# spdlog::spdlog
# PRIVATE
# igl::cgal
# )

# set_target_properties(libTetWild PROPERTIES OUTPUT_NAME "tetwild")
# target_link_libraries(TetWild
# libTetWild
# CLI11::CLI11
# igl::cgal
# )

# ################################################################################
# # Folders for Visual Studio/XCode IDEs
# ################################################################################

# # geogram
# set_target_properties(geogram PROPERTIES FOLDER extern/geogram)
# set_target_properties(geogram_third_party PROPERTIES FOLDER extern/geogram)
# set_target_properties(uninstall PROPERTIES FOLDER extern/geogram)

# # pymesh
# set_target_properties(pymesh_tiny PROPERTIES FOLDER extern/pymesh)