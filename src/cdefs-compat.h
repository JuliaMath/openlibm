#ifndef _CDEFS_COMPAT_H_
#define	_CDEFS_COMPAT_H_

#if defined(__cplusplus)
#define	__BEGIN_DECLS	extern "C" {
#define	__END_DECLS	}
#else
#define	__BEGIN_DECLS
#define	__END_DECLS
#endif

#endif /* _CDEFS_COMPAT_H_ */
