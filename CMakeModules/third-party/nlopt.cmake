
if(TARGET NLopt::nlopt)
  return()
endif()

message(STATUS "Loading Nlopt...")

set(NLOPT_CXX ON CACHE INTERNAL "enable cxx routines" FORCE)
set(NLOPT_FORTRAN OFF CACHE INTERNAL "enable fortran" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "Build NLopt as a shared library" FORCE)
set(NLOPT_PYTHON OFF CACHE INTERNAL "build python bindings" FORCE)
set(NLOPT_OCTAVE OFF CACHE INTERNAL "build octave bindings" FORCE)
set(NLOPT_MATLAB OFF CACHE INTERNAL "build matlab bindings" FORCE)
set(NLOPT_GUILE OFF CACHE INTERNAL "build guile bindings" FORCE)
set(NLOPT_JAVA OFF CACHE INTERNAL "build java bindings" FORCE)
set(NLOPT_SWIG OFF CACHE INTERNAL "use SWIG to build bindings" FORCE)
set(NLOPT_LUKSAN ON CACHE INTERNAL "enable LGPL Luksan solvers" FORCE)
set(NLOPT_TESTS OFF CACHE INTERNAL "build unit tests" FORCE)


include(FetchContent)
FetchContent_Declare(
  NLopt
  URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.10.0.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS NAMES NLopt
)

FetchContent_MakeAvailable(NLopt)

message(STATUS "Done.")
