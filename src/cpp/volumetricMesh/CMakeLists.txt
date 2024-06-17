set(VOLUMETRIC_MESH_H
  cubicMesh.h
  tetMesh.h
  volumetricMeshENuMaterial.h
  volumetricMesh.h
  volumetricMeshMooneyRivlinMaterial.h
  volumetricMeshOrthotropicMaterial.h
  volumetricMeshParser.h
  generateMassMatrix.h
)
set(VOLUMETRIC_MESH_SRC
  cubicMesh.cpp
  tetMesh.cpp
  volumetricMesh.cpp
  volumetricMeshENuMaterial.cpp
  volumetricMeshMooneyRivlinMaterial.cpp
  volumetricMeshOrthotropicMaterial.cpp
  volumetricMeshParser.cpp
  generateMassMatrix.cpp
)

if(NOT TARGET mesh)
  message(STATUS "volumetricMesh is not included. Missing: mesh")
  return()
endif()

add_libpgo_lib(volumetricMesh "${VOLUMETRIC_MESH_SRC}" "${VOLUMETRIC_MESH_H}")

target_link_libraries(volumetricMesh PUBLIC mesh)