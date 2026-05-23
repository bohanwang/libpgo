if(TARGET pgoLogging)
  add_libpgo_tools(computeStressFieldStats computeStressFieldStats.cpp "pgoLogging")
  target_link_libraries(computeStressFieldStats PRIVATE argparse::argparse nlohmann_json::nlohmann_json fmt::fmt)
endif()
