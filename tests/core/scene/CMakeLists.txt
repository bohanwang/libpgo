set(AVAILABLE_LIBS)

foreach(lib ${PGO_GLOBAL_LIBRARY_TARGETS})
  if(TARGET ${lib})
    list(APPEND AVAILABLE_LIBS ${lib})
  endif()
endforeach()

add_libpgo_tools(core_scene_test core_scene_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_scene_test PRIVATE GTest::gtest_main)
set_property(TARGET core_scene_test PROPERTY FOLDER "tests")

add_test(NAME core.scene.geometry COMMAND core_scene_test)

add_libpgo_tools(core_scene_contact_embedding_test contact_embedding_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_scene_contact_embedding_test PRIVATE GTest::gtest_main)
set_property(TARGET core_scene_contact_embedding_test PROPERTY FOLDER "tests")

add_test(NAME core.scene.contact_embedding COMMAND core_scene_contact_embedding_test)

add_libpgo_tools(core_scene_self_contact_detection_test self_contact_detection_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_scene_self_contact_detection_test PRIVATE GTest::gtest_main)
set_property(TARGET core_scene_self_contact_detection_test PROPERTY FOLDER "tests")

add_test(NAME core.scene.self_contact_detection COMMAND core_scene_self_contact_detection_test)

add_libpgo_tools(core_scene_self_contact_handler_test self_contact_handler_test.cpp "${AVAILABLE_LIBS}")
target_link_libraries(core_scene_self_contact_handler_test PRIVATE GTest::gtest_main)
set_property(TARGET core_scene_self_contact_handler_test PROPERTY FOLDER "tests")

add_test(NAME core.scene.self_contact_handler COMMAND core_scene_self_contact_handler_test)
