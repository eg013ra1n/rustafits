#ifndef COMPAT_H
#define COMPAT_H

// Thread-local storage
#ifdef _MSC_VER
  #define THREAD_LOCAL __declspec(thread)
#else
  #define THREAD_LOCAL __thread
#endif

// Byte-swap intrinsics
#ifdef _MSC_VER
  #include <stdlib.h>
  #define BSWAP16(x) _byteswap_ushort(x)
  #define BSWAP32(x) _byteswap_ulong(x)
  #define BSWAP64(x) _byteswap_uint64(x)
#else
  #define BSWAP16(x) __builtin_bswap16(x)
  #define BSWAP32(x) __builtin_bswap32(x)
  #define BSWAP64(x) __builtin_bswap64(x)
#endif

// Case-insensitive string compare
#ifdef _MSC_VER
  #define STRCASECMP _stricmp
#else
  #define STRCASECMP strcasecmp
#endif

#endif
