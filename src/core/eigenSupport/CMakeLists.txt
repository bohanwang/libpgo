set(EIGENSUPPORT_HEADERS
  EigenSupport.h
  EigenDef.h
)

set(EIGENSUPPORT_SOURCES
  EigenSupport.cpp
)

if(TARGET MKL::MKL)
  list(APPEND EIGENSUPPORT_HEADERS EigenMKLPardisoSupport.h)
  list(APPEND EIGENSUPPORT_SOURCES EigenMKLPardisoSupport.cpp)
endif()

if(PGO_HAS_ARPACK)
  message(STATUS "has arpack")
  list(APPEND EIGENSUPPORT_HEADERS EigenARPACKSupport.h)
  list(APPEND EIGENSUPPORT_SOURCES EigenARPACKSupport.cpp)
else()
  message(STATUS "has no arpack")
endif()

set(EIGENSUPPORT_DEPS Eigen3::Eigen TBB::tbb nlohmann_json::nlohmann_json)

foreach(tgt ${EIGENSUPPORT_DEPS})
  if(NOT TARGET ${tgt})
    message(STATUS "eigenSupport is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(eigenSupport "${EIGENSUPPORT_SOURCES}" "${EIGENSUPPORT_HEADERS}")

target_link_libraries(eigenSupport PUBLIC ${EIGENSUPPORT_DEPS})

if(PGO_HAS_ARPACK)
  target_link_libraries(eigenSupport PUBLIC arpack ${GFORTRAN_LIBRARY} ${QUADMATH_LIBRARY})
  target_compile_definitions(eigenSupport PUBLIC PGO_HAS_ARPACK)
endif()
