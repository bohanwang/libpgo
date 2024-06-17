set(MESH_HEADERS triMeshGeo.h createTriMesh.h triMeshNeighbor.h triMeshPseudoNormal.h tetMeshGeo.h edgeKey.h triKey.h rectKey.h tetKey.h
  triMeshManifold.h tetMeshManifold.h tetrahedron.h windingNumberTree.h labelOuterTets.h geometryQuery.h boundingBox.h simpleSphere.h
  triangle.h tribox3.h halfSpace.h predicates.h initPredicates.h verticesInfo.h meshIntersection.h plane.h boundingVolumeTreeBase.h
  boundingVolumeTree.h verticesCorrespondence.h pointInsideOutsideQuery.h triangleSampler.h voxelFunction.h triMeshSampler.h meshLinearAlgebra.h)

set(MESH_SOURCES triMeshGeo.cpp createTriMesh.cpp triMeshNeighbor.cpp triMeshPseudoNormal.cpp tetMeshGeo.cpp triKey.cpp rectKey.cpp tetKey.cpp
  triMeshManifold.cpp tetMeshManifold.cpp tetrahedron.cpp windingNumberTree.cpp labelOuterTets.cpp geometryQuery.cpp boundingBox.cpp simpleSphere.cpp
  triangle.cpp tribox3.cpp predicates.cpp verticesInfo.cpp meshIntersection.cpp boundingVolumeTreeBase.cpp boundingVolumeTree.cpp
  verticesCorrespondence.cpp pointInsideOutsideQuery.cpp lineSegmentDistance.cpp triMeshSampler.cpp triMeshGeoLoader.cpp)

set(MESH_DEPS basicAlgorithms eigenSupport macros pgoLogging)

foreach(tgt ${MESH_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "mesh is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

if(NOT TARGET CGAL::CGAL)
  list(APPEND MESH_SOURCES intersection_tunicate.cpp predicates_Shewchuk.c)
endif()

add_libpgo_lib(mesh "${MESH_SOURCES}" "${MESH_HEADERS}")

target_link_libraries(mesh PUBLIC ${MESH_DEPS})

# use cgal predicates, the old code has some problems when running on modern hardware
# intersection_tunicate.cpp predicates_Shewchuk.c
if(TARGET CGAL::CGAL)
  target_compile_definitions(mesh PUBLIC PGO_USE_CGAL_PREDICATES)
  target_link_libraries(mesh PUBLIC CGAL::CGAL)
endif()

target_link_libraries(mesh PUBLIC TBB::tbb)
target_link_libraries(mesh PUBLIC tiny_obj_loader)
