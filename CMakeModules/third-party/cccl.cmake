if(TARGET CCCL::CCCL)
  return()
endif()

message(STATUS "Loading cccl...")

include(FetchContent)
FetchContent_Declare(
  cccl
  GIT_REPOSITORY https://github.com/NVIDIA/cccl.git
  GIT_TAG 8306d992387dfafa9e7d81387b4bf2c1c30ee8bf
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS 3.1.0
)

FetchContent_MakeAvailable(cccl)

message(STATUS "Done.")


