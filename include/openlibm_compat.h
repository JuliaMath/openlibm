#ifndef _OPENLIBM_COMPAT_H_
#define	_OPENLIBM_COMPAT_H_

#define DLLEXPORT

#define __FBSDID(s)

#ifdef __MINIOS__
/* No stdio.h on Mini-OS. */
#include <sys/cdefs.h>
#else
/*
 * We cannot be certain that this operating system has <sys/cdefs.h>.
 * Instead, include a header file that is likely to pull in this header.
 */
#include <stdio.h>
#endif

#if defined(__cplusplus)
#define	__BEGIN_DECLS	extern "C" {
#define	__END_DECLS	}
#else
#define	__BEGIN_DECLS
#define	__END_DECLS
#endif

#ifdef __GNUC__
#ifndef __strong_reference
#ifdef __APPLE__
#define __strong_reference(sym,aliassym) __weak_reference(sym,aliassym)
#else
#define __strong_reference(sym,aliassym)	\
	DLLEXPORT extern __typeof (sym) aliassym __attribute__ ((__alias__ (#sym)));
#endif /* __APPLE__ */
#endif /* __strong_reference */

#ifndef __weak_reference
#ifdef __ELF__
#ifdef __STDC__
#define	__weak_reference(sym,alias)	\
	__asm__(".weak " #alias);	\
	__asm__(".equ "  #alias ", " #sym)
#define	__warn_references(sym,msg)	\
	__asm__(".section .gnu.warning." #sym);	\
	__asm__(".asciz \"" msg "\"");	\
	__asm__(".previous")
#else
#define	__weak_reference(sym,alias)	\
	__asm__(".weak alias");		\
	__asm__(".equ alias, sym")
#define	__warn_references(sym,msg)	\
	__asm__(".section .gnu.warning.sym"); \
	__asm__(".asciz \"msg\"");	\
	__asm__(".previous")
#endif	/* __STDC__ */
#elif defined(__clang__) /* CLANG */
#ifdef __STDC__
#define __weak_reference(sym,alias)     \
    __asm__(".weak_reference " #alias); \
    __asm__(".set " #alias ", " #sym)
#else
#define __weak_reference(sym,alias)     \
    __asm__(".weak_reference alias");\
    __asm__(".set alias, sym")
#endif
#else	/* !__ELF__ */
#ifdef __STDC__
#define __weak_reference(sym,alias)	\
	__asm__(".stabs \"_" #alias "\",11,0,0,0");	\
	__asm__(".stabs \"_" #sym "\",1,0,0,0")
#define __warn_references(sym,msg)	\
	__asm__(".stabs \"" msg "\",30,0,0,0");		\
	__asm__(".stabs \"_" #sym "\",1,0,0,0")
#else
#define __weak_reference(sym,alias)	\
	__asm__(".stabs \"_/**/alias\",11,0,0,0");	\
	__asm__(".stabs \"_/**/sym\",1,0,0,0")
#define __warn_references(sym,msg)	\
	__asm__(".stabs msg,30,0,0,0");			\
	__asm__(".stabs \"_/**/sym\",1,0,0,0")
#endif	/* __STDC__ */
#endif	/* __ELF__ */
#endif  /* __weak_reference */
#endif	/* __GNUC__ */


#include <stdint.h>
#include <limits.h>
#include <stdint.h>

#ifdef __GLIBC__
/* Not sure what to do about __pure2 on linux */
#define __pure2
#endif

#ifdef _WIN32
/* Not sure what to do about __pure2 on windows */
#define __pure2
#endif

typedef uint8_t               u_int8_t;
typedef uint16_t              u_int16_t;
typedef uint32_t              u_int32_t;
typedef uint64_t              u_int64_t;

typedef float float_t;

#ifndef __always_inline
# define __always_inline __inline __attribute__ ((__always_inline__))
#endif

#endif
