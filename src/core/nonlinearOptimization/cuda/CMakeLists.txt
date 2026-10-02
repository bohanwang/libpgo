set(CUDA_OPTIMIZATION_HEADERS
  cudaPotentialEnergy.h
  cudaBlas.h
  cudaNewtonSolver.h
  cudaFrozenHessianEnergy.h
  cudaQuadraticEnergy.h)
set(CUDA_OPTIMIZATION_SOURCES
  cudaNewtonSolver.cu
  cudaFrozenHessianEnergy.cu
  cudaQuadraticEnergy.cu)

add_libpgo_lib(cudaNonlinearOptimization "${CUDA_OPTIMIZATION_SOURCES}" "${CUDA_OPTIMIZATION_HEADERS}")
target_link_libraries(cudaNonlinearOptimization PUBLIC cudaUtilities eigenSupport PRIVATE CUDA::cublas)
target_include_directories(cudaNonlinearOptimization PRIVATE "${PROJECT_SOURCE_DIR}/third-party/libuipc/include")
target_compile_options(cudaNonlinearOptimization PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
set_target_properties(cudaNonlinearOptimization PROPERTIES CUDA_STANDARD_REQUIRED ON)
