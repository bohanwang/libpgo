add_libpgo_lib(cudaContact
  "cudaEmbeddedSurfaceIPCPotentialEnergy.cu;cudaSurfaceIPCBroadPhase.cu;cudaSurfaceIPCFrictionPotentialEnergy.cu"
  "cudaEmbeddedSurfaceIPCPotentialEnergy.h;cudaSurfaceIPCBroadPhase.h;cudaSurfaceIPCFrictionPotentialEnergy.h;cudaIPCGeometry.h;cudaIPCContinuousCollision.h")
target_link_libraries(cudaContact PUBLIC cudaNonlinearOptimization contact)
target_include_directories(cudaContact PRIVATE "${PROJECT_SOURCE_DIR}/third-party/libuipc/include")
# CPU BLAS specializations do not implement device-side Eigen eigensolves.
# Keep the imported libuipc 12x12 PSD projection on Eigen's CUDA path.
target_compile_options(cudaContact PRIVATE
  "$<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr;-UEIGEN_USE_MKL_ALL;-UEIGEN_USE_BLAS;-UEIGEN_USE_LAPACKE>")
set_target_properties(cudaContact PROPERTIES CUDA_STANDARD_REQUIRED ON)
