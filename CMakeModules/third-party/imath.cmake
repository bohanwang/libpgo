if(TARGET Imath::Imath)
  return()
endif()

message(STATUS "Loading Imath...")

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build third-party libraries statically" FORCE)
set(IMATH_INSTALL ON CACHE BOOL "Define Imath install/export targets" FORCE)
set(IMATH_INSTALL_PKG_CONFIG OFF CACHE BOOL "Install Imath pkg-config file" FORCE)
set(PYTHON OFF CACHE BOOL "Build Imath Boost.Python bindings" FORCE)
set(PYBIND11 OFF CACHE BOOL "Build Imath pybind11 bindings" FORCE)

include(FetchContent)
FetchContent_Declare(
  Imath
  URL https://github.com/AcademySoftwareFoundation/Imath/archive/refs/tags/v3.2.2.tar.gz
  URL_HASH SHA256=b4275d83fb95521510e389b8d13af10298ed5bed1c8e13efd961d91b1105e462
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  OVERRIDE_FIND_PACKAGE
)

FetchContent_MakeAvailable(Imath)

# Imath 3.2.2 exposes src/Imath as its build-tree include directory, while
# consumers (including Alembic) use installed-style includes such as
# <Imath/half.h>. Add the parent directory for build-tree consumers only.
target_include_directories(
  Imath
  INTERFACE $<BUILD_INTERFACE:${imath_SOURCE_DIR}/src>
)

message(STATUS "Done.")
