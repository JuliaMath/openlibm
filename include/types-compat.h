#ifndef _TYPES_COMPAT_H_
#define	_TYPES_COMPAT_H_

#if defined(_MSC_VER) && !defined(_WIN32)
    #define _WIN32
#endif

#include <sys/types.h>
#include <limits.h>
#include <stdint.h>

#ifdef __linux__
/* Not sure what to do about __pure2 on linux */
#define __pure2 
#endif

#ifdef _WIN32
/* Not sure what to do about __pure2 on windows */
#define __pure2 
typedef uint8_t               u_int8_t;
typedef uint16_t              u_int16_t;
typedef uint32_t              u_int32_t;
typedef uint64_t              u_int64_t;
#endif


#endif
