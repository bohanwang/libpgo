add_libpgo_lib(cudaSolidDeformationModel
  "cudaDeformationModelEnergy.cu"
  "cudaDeformationModelEnergy.h;cudaHexJacobian.h")
add_libpgo_lib(cudaKoiterDeformationModel
  "cudaKoiterDeformationModelEnergy.cu"
  "cudaKoiterDeformationModelEnergy.h")
# CPU Eigen links MKL through its imported target. Device constitutive algebra
# must use Eigen's CUDA implementation; NVCC 13.1 also crashes when optimizing
# projected tangents with the inherited CPU BLAS specializations.
foreach(target cudaSolidDeformationModel cudaKoiterDeformationModel)
  target_link_libraries(${target} PUBLIC cudaNonlinearOptimization solidDeformationModel)
  target_include_directories(${target} PRIVATE "${PROJECT_SOURCE_DIR}/third-party/libuipc/include")
  target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
  target_compile_options(${target} PRIVATE
    "$<$<COMPILE_LANGUAGE:CUDA>:-UEIGEN_USE_MKL_ALL;-UEIGEN_USE_BLAS;-UEIGEN_USE_LAPACKE>")
  set_target_properties(${target} PROPERTIES CUDA_STANDARD_REQUIRED ON)
endforeach()
