file(GLOB PROJECTS_ALL LIST_DIRECTORIES TRUE *)

foreach(prj ${PROJECTS_ALL})
  message(STATUS "checking ${prj} ...")

  if(IS_DIRECTORY ${prj})
    add_subdirectory(${prj})
  endif()
endforeach(prj ${PROJECTS_ALL})
