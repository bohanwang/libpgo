#pragma once

#ifndef _LIBPGO_C_DEF_H_
#  define _LIBPGO_C_DEF_H_

#  ifdef __cplusplus
#    define LIBPGO_EXTERN_C_BEG \
      extern "C"                \
      {
#    define LIBPGO_EXTERN_C_END }
#  else
#    define LIBPGO_EXTERN_C_BEG
#    define LIBPGO_EXTERN_C_END
#  endif

#  ifdef LIBPGO_C_BUILD_DLL
#    ifdef _MSC_VER
#      define LIBPGO_C_EXPORT __declspec(dllexport)
#    else
#      define LIBPGO_C_EXPORT __attribute__((visibility("default")))
#    endif
#  elif LIBPGO_C_LOAD_DLL
#    ifdef _MSC_VER
#      define LIBPGO_C_EXPORT __declspec(dllimport)
#    else
#      define LIBPGO_C_EXPORT
#    endif
#  else
#    define LIBPGO_C_EXPORT
#  endif
#endif