#ifndef _TYPES_COMPAT_H_
#define	_TYPES_COMPAT_H_

#if (defined(_WIN32) || defined (_MSC_VER)) && !defined(__WIN32__)
    #define __WIN32__
#endif

#ifdef __APPLE__
#include <sys/_types.h>
#include <machine/_limits.h>
#endif

#ifdef __linux__
/* Not sure what to do about __pure2 on linux */
#define __pure2 
#include <sys/types.h>
#include <limits.h>
#endif

#ifdef __WIN32__
/* Not sure what to do about __pure2 on linux */
#define __pure2 
#include <sys/types.h>
#include <stdint.h>
typedef uint8_t               u_int8_t;
typedef uint16_t              u_int16_t;
typedef uint32_t              u_int32_t;
typedef uint64_t              u_int64_t;
#include <limits.h>
#endif


#endif