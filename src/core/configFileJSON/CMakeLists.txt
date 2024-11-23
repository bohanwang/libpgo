set(CONFIGFILEJSON_DEPS pgoLogging nlohmann_json::nlohmann_json)

foreach (tgt ${CONFIGFILEJSON_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "configFileJSON is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(configFileJSON configFileJSON.cpp configFileJSON.h)

target_link_libraries(configFileJSON PUBLIC ${CONFIGFILEJSON_DEPS})