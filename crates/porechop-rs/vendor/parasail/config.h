#ifndef PORECHOP_PARASAIL_CONFIG_H
#define PORECHOP_PARASAIL_CONFIG_H

#include <limits.h>
#if INT_MAX == 2147483647
#define SIZEOF_INT 4
#else
#error "Parasail requires 32-bit int"
#endif

#ifdef _WIN32
#define HAVE__ALIGNED_MALLOC 1
#include <malloc.h>
#else
#define HAVE_POSIX_MEMALIGN 1
#endif

#endif
