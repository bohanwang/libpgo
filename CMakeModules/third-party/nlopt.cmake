
if(TARGET NLopt::nlopt)
  return()
endif()

message(STATUS "Loading Nlopt...")

pgo_dep_option(NLOPT_CXX INTERNAL ON "enable cxx routines")
pgo_dep_option(NLOPT_FORTRAN INTERNAL OFF "enable fortran")
pgo_dep_option(BUILD_SHARED_LIBS INTERNAL OFF "Build NLopt as a shared library")
pgo_dep_option(NLOPT_PYTHON INTERNAL OFF "build python bindings")
pgo_dep_option(NLOPT_OCTAVE INTERNAL OFF "build octave bindings")
pgo_dep_option(NLOPT_MATLAB INTERNAL OFF "build matlab bindings")
pgo_dep_option(NLOPT_GUILE INTERNAL OFF "build guile bindings")
pgo_dep_option(NLOPT_JAVA INTERNAL OFF "build java bindings")
pgo_dep_option(NLOPT_SWIG INTERNAL OFF "use SWIG to build bindings")
pgo_dep_option(NLOPT_LUKSAN INTERNAL ON "enable LGPL Luksan solvers")
pgo_dep_option(NLOPT_TESTS INTERNAL OFF "build unit tests")


include(FetchContent)
FetchContent_Declare(
  NLopt
  URL https://github.com/stevengj/nlopt/archive/refs/tags/v2.10.0.zip
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

pgo_fetch_make_available(NLopt)

message(STATUS "Done.")
