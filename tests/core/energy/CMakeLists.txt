set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_libpgo_tools(core_energy_material_binding_test simulationMesh_material_binding_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_energy_material_binding_test PRIVATE GTest::gtest_main)
set_property(TARGET core_energy_material_binding_test PROPERTY FOLDER "tests")

add_test(NAME core.energy.material_binding COMMAND core_energy_material_binding_test)

add_libpgo_tools(core_energy_cubic_deformation_model_test cubicMeshDeformationModel_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_energy_cubic_deformation_model_test PRIVATE GTest::gtest_main)
set_property(TARGET core_energy_cubic_deformation_model_test PROPERTY FOLDER "tests")

add_test(NAME core.energy.cubic_deformation_model COMMAND core_energy_cubic_deformation_model_test)

add_libpgo_tools(core_energy_point_penetration_barrier_test pointPenetrationBarrierEnergy_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_energy_point_penetration_barrier_test PRIVATE GTest::gtest_main)
set_property(TARGET core_energy_point_penetration_barrier_test PROPERTY FOLDER "tests")

add_test(NAME core.energy.point_penetration_barrier COMMAND core_energy_point_penetration_barrier_test)

add_libpgo_tools(core_energy_point_triangle_pair_barrier_test pointTrianglePairBarrierEnergy_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_energy_point_triangle_pair_barrier_test PRIVATE GTest::gtest_main)
set_property(TARGET core_energy_point_triangle_pair_barrier_test PROPERTY FOLDER "tests")

add_test(NAME core.energy.point_triangle_pair_barrier COMMAND core_energy_point_triangle_pair_barrier_test)

add_libpgo_tools(core_energy_newton_raphson_solver_test newtonRaphsonSolver_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_energy_newton_raphson_solver_test PRIVATE GTest::gtest_main)
set_property(TARGET core_energy_newton_raphson_solver_test PROPERTY FOLDER "tests")

add_test(NAME core.energy.newton_raphson_solver COMMAND core_energy_newton_raphson_solver_test)
