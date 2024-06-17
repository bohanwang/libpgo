set(SOLID_DEFORMATION_MODEL_HEADERS 
  elasticModel.h
  elasticModelInvariantBasedMaterial.h
  elasticModelStableNeoHookeanMaterial.h
  elasticModelLinearMaterial.h
  elasticModelHillTypeMaterial.h
  elasticModelVolumeMaterial.h
  elasticModelCombinedMaterial.h

  elasticModel1DQuadratic.h
  elasticModel1DCubicSpline.h

  invariantBasedMaterial.h
  invariantBasedMaterialStVK.h
  invariantBasedMaterialNeoHookean.h

  plasticModel.h
  plasticModel3DConstant.h
  plasticModel3D3DOF.h
  plasticModel3D6DOF.h

  deformationModel.h
  tetMeshDeformationModel.h
  deformationModelManager.h

  deformationModelEnergy.h
  deformationModelAssembler.h
  deformationModelAssemblerImpl.h
  deformationModelAssemblerElastic.h
  deformationModelAssemblerElasticImpl.h
  deformationModelAssemblerElasticFullPlastic.h
  deformationModelAssemblerElasticFullPlasticImpl.h
  deformationModelAssemblerElasticMultiPlastic.h
  deformationModelAssemblerElasticMultiPlasticImpl.h
  hessianMatrixHandle.h

  simulationMesh.h
  #deformationModelFDTest.h

  tetVolumeConstraintFunctions.h
  tetMeshOccupation.h

  segmentChainConstraintFunctions.h
  segmentBinormalConstraintFunctions.h

  prescribedPrincipleStressConstraintFunctions.h

  generateTetMeshMatrix.h
)

set(SOLID_DEFORMATION_MODEL_SOURCES 
  elasticModelInvariantBasedMaterial.cpp
  elasticModelStableNeoHookeanMaterial.cpp
  elasticModelLinearMaterial.cpp
  elasticModelHillTypeMaterial.cpp
  elasticModelVolumeMaterial.cpp

  elasticModel1DQuadratic.cpp
  elasticModel1DCubicSpline.cpp

  invariantBasedMaterialStVK.cpp
  invariantBasedMaterialNeoHookean.cpp
  
  plasticModel.cpp
  plasticModel3DConstant.cpp
  plasticModel3D3DOF.cpp
  plasticModel3D6DOF.cpp

  tetMeshDeformationModel.cpp
  deformationModelManager.cpp

  deformationModelEnergy.cpp
  deformationModelAssembler.cpp
  deformationModelAssemblerElastic.cpp
  deformationModelAssemblerElasticMultiPlastic.cpp
  deformationModelAssemblerElasticFullPlastic.cpp

  simulationMesh.cpp
  #deformationModelFDTest.cpp
  
  tetVolumeConstraintFunctions.cpp
  tetMeshOccupation.cpp

  segmentChainConstraintFunctions.cpp
  segmentBinormalConstraintFunctions.cpp

  prescribedPrincipleStressConstraintFunctions.cpp

  generateTetMeshMatrix.cpp
)

set(SOLID_DEFORMATION_MODEL_DEPS 
  pgoLogging eigenSupport nonlinearOptimization configFileJSON volumetricMesh 
  Eigen3::Eigen TBB::tbb fmt::fmt-header-only 
)

foreach(tgt ${SOLID_DEFORMATION_MODEL_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "solidDeformationModel is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(solidDeformationModel "${SOLID_DEFORMATION_MODEL_SOURCES}" "${SOLID_DEFORMATION_MODEL_HEADERS}")

target_link_libraries(solidDeformationModel PUBLIC ${SOLID_DEFORMATION_MODEL_DEPS})