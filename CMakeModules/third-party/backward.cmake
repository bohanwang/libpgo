include(FetchContent)

# Also requires one of: libbfd (gnu binutils), libdwarf, libdw (elfutils)
FetchContent_Declare(backward
    GIT_REPOSITORY https://github.com/bombela/backward-cpp
    GIT_TAG 0bfd0a07a61551413ccd2ab9a9099af3bad40681
    SYSTEM          # optional, the Backward include directory will be treated as system directory
)
FetchContent_MakeAvailable(backward)
