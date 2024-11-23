set(CGALINTERFACE_HEADERS cgalBasic.h cgalInterface.h cgalTemplateUtilities.h)
set(CGALINTERFACE_SOURCES cgalInterface.cpp)
set(CGALINTERFACE_DEPS pgoLogging eigenSupport mesh CGAL::CGAL CGAL::TBB_support)

foreach (tgt ${CGALINTERFACE_DEPS})
  if (NOT TARGET ${tgt})
    message(STATUS "cgalInterface is not included. Missing: ${tgt}")
    return()
  endif()
endforeach()

add_libpgo_lib(cgalInterface "${CGALINTERFACE_SOURCES}" "${CGALINTERFACE_HEADERS}")

target_link_libraries(cgalInterface PUBLIC ${CGALINTERFACE_DEPS})

if (TARGET Ceres::ceres)
  target_link_libraries(cgalInterface PUBLIC Ceres::ceres)
endif()