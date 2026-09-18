if(TARGET autodiff::autodiff)
  return()
endif()

message(STATUS "Loading autodiff...")

set(AUTODIFF_BUILD_TESTS OFF CACHE BOOL "Enable the compilation of the test files." FORCE)
set(AUTODIFF_BUILD_PYTHON OFF CACHE BOOL "Enable the compilation of the python bindings." FORCE)
set(AUTODIFF_BUILD_EXAMPLES OFF CACHE BOOL "Enable the compilation of the example files." FORCE)
set(AUTODIFF_BUILD_DOCS OFF CACHE BOOL "Enable the build of the documentation and website." FORCE)

include(FetchContent)
FetchContent_Declare(
  autodiff_pkg
  URL https://github.com/autodiff/autodiff/archive/refs/tags/v1.1.2.tar.gz
  URL_HASH SHA256=86f68aabdae1eed214bfbf0ddaa182c78ea1bb99e4df404efb7b94d30e06b744
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
)

FetchContent_MakeAvailable(autodiff_pkg)

message(STATUS "Done.")
