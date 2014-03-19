#ifndef _CDEFS_COMPAT_H_
#define	_CDEFS_COMPAT_H_

#if defined(_MSC_VER) && !defined(_WIN32)
    #define _WIN32
#endif

#ifndef _WIN32
#include "sys/cdefs.h"
#else /* _WIN32 */

#if defined(__cplusplus)
#define	__BEGIN_DECLS	extern "C" {
#define	__END_DECLS	}
#else
#define	__BEGIN_DECLS
#define	__END_DECLS
#endif

#define _SYS_CDEFS_H_

#endif /* _WIN32 */




#ifdef __GNUC__
#ifndef __strong_reference
#define __strong_reference(sym,aliassym)	
	//extern __typeof (sym) aliassym __attribute__ ((__alias__ (#sym)));
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


#endif /* _CDEFS_COMPAT_H_ */
