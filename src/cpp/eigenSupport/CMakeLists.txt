set(EIGENSUPPORT_HEADERS 
  EigenSupport.h
  EigenDef.h
)

set(EIGENSUPPORT_SOURCES 
  EigenSupport.cpp
)

if (TARGET ARPACK::ARPACK)
  list(APPEND EIGENSUPPORT_HEADERS EigenARPACKSupport.h)
  list(APPEND EIGENSUPPORT_SOURCES EigenARPACKSupport.cpp)
endif()

if (TARGET MKL::MKL)
  list(APPEND EIGENSUPPORT_HEADERS EigenMKLPardisoSupport.h)
  list(APPEND EIGENSUPPORT_SOURCES EigenMKLPardisoSupport.cpp)
endif()

set(EIGENSUPPORT_DEPS Eigen3::Eigen TBB::tbb nlohmann_json::nlohmann_json)

foreach (tgt ${EIGENSUPPORT_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "eigenSupport is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(eigenSupport "${EIGENSUPPORT_SOURCES}" "${EIGENSUPPORT_HEADERS}")

target_link_libraries(eigenSupport PUBLIC ${EIGENSUPPORT_DEPS})