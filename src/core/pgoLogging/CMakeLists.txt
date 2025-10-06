set(PGOLOGGING_HEADERS pgoLogging.h fmtEigen.h)
set(PGOLOGGING_DEPS Eigen3::Eigen spdlog::spdlog_header_only fmt::fmt-header-only)

foreach(tgt ${PGOLOGGING_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "pgoLogging is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(pgoLogging pgoLogging.cpp "${PGOLOGGING_HEADERS}")

target_link_libraries(pgoLogging PUBLIC ${PGOLOGGING_DEPS})
target_compile_definitions(pgoLogging PUBLIC FMT_UNICODE=0)
