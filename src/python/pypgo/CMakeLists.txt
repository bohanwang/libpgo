find_package(Python COMPONENTS Interpreter Development REQUIRED)

nanobind_add_module(
  pypgo
  NB_DOMAIN pypgo
  pypgo.cpp
)

target_link_libraries(pypgo PRIVATE
  nonlinearOptimization
  contact
  solidDeformationModel
  simulation
  mesh
)

target_compile_definitions(pypgo PRIVATE NB_EIGEN)
