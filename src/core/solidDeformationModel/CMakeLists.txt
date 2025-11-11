set(SOLID_DEFORMATION_MODEL_HEADERS
  elasticModel.h
  elasticModel3DDeformationGradient.h
  elasticModelStableNeoHookeanMaterial.h
  elasticModelInvariantBasedMaterial.h  
  elasticModelLinearMaterial.h
  elasticModelHillTypeMaterial.h
  elasticModelVolumeMaterial.h
  elasticModelCombinedMaterial.h
  elasticModel3DMooneyRivlin.h
  elasticModel3DSTVKMaterial.h

  # elasticModel1DQuadratic.h
  # elasticModel1DCubicSpline.h

  elasticModel2DFundamentalForms.h
  elasticModel2DFundamentalFormsFabric.h
  elasticModel2DFundamentalFormsSTVK.h

  invariantBasedMaterial.h
  invariantBasedMaterialStVK.h
  invariantBasedMaterialNeoHookean.h

  plasticModel.h
  plasticModel3DDeformationGradient.h
  plasticModel3DConstant.h
  plasticModel3D3DOF.h
  plasticModel3D6DOF.h

  plasticModel2DFundamentalForms.h
  plasticModel2DFundamentalFormsUniformStretch.h
  plasticModel2DFundamentalFormsFirstOnly.h

  deformationModel.h
  tetMeshDeformationModel.h  
  koiterDeformationModel.h

  deformationModelManager.h

  deformationModelEnergy.h
  deformationModelAssembler.h
  hessianMatrixHandle.h

  simulationMesh.h
  deformationModelFDTest.h

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
  elasticModel3DMooneyRivlin.cpp
  elasticModel3DSTVKMaterial.cpp

  # elasticModel1DQuadratic.cpp
  # elasticModel1DCubicSpline.cpp  
  
  elasticModel2DFundamentalFormsFabric.cpp
  elasticModel2DFundamentalFormsSTVK.cpp

  invariantBasedMaterialStVK.cpp
  invariantBasedMaterialNeoHookean.cpp
  
  plasticModel3DConstant.cpp
  plasticModel3D3DOF.cpp
  plasticModel3D6DOF.cpp

  plasticModel2DFundamentalForms.cpp
  plasticModel2DFundamentalFormsUniformStretch.cpp
  plasticModel2DFundamentalFormsFirstOnly.cpp

  tetMeshDeformationModel.cpp
  koiterDeformationModel.cpp

  deformationModelManager.cpp

  deformationModelEnergy.cpp
  deformationModelAssembler.cpp

  simulationMesh.cpp
  deformationModelFDTest.cpp

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
