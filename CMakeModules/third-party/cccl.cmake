if(TARGET CCCL::CCCL)
  return()
endif()

message(STATUS "Loading cccl...")

include(FetchContent)
FetchContent_Declare(
  cccl
  # GIT_REPOSITORY https://github.com/NVIDIA/cccl.git
  # GIT_TAG 8306d992387dfafa9e7d81387b4bf2c1c30ee8bf
  URL https://github.com/NVIDIA/cccl/releases/download/v3.1.3/cccl-src-v3.1.3.zip
  URL_HASH SHA256=ca8cbb65bb9f0d8d9734d5f70a193f243e088cb5f4594c7d40b09fae34ae3c18
  EXCLUDE_FROM_ALL
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  FIND_PACKAGE_ARGS 3.1.3
)

FetchContent_MakeAvailable(cccl)

message(STATUS "Done.")

