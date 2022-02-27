/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

/*
 * from: @(#)fdlibm.h 5.1 93/09/24
 * $FreeBSD: src/lib/msun/src/math_private.h,v 1.34 2011/10/21 06:27:56 das Exp $
 */

#ifndef _MATH_PRIVATE_H_
#define	_MATH_PRIVATE_H_

#include <openlibm_complex.h>
#include <openlibm_defs.h>
#include "cdefs-compat.h"
#include "types-compat.h"
#include "fpmath.h"
#include <stdint.h>
#include "math_private_openbsd.h"

/*
 * The original fdlibm code used statements like:
 *	n0 = ((*(int*)&one)>>29)^1;		* index of high word *
 *	ix0 = *(n0+(int*)&x);			* high word of x *
 *	ix1 = *((1-n0)+(int*)&x);		* low word of x *
 * to dig two 32 bit words out of the 64 bit IEEE floating point
 * value.  That is non-ANSI, and, moreover, the gcc instruction
 * scheduler gets it wrong.  We instead use the following macros.
 * Unlike the original code, we determine the endianness at compile
 * time, not at run time; I don't see much benefit to selecting
 * endianness at run time.
 */

/*
 * A union which permits us to convert between a double and two 32 bit
 * ints.
 */

#if __FLOAT_WORD_ORDER__ == __ORDER_BIG_ENDIAN__

typedef union
{
  double value;
  struct
  {
    u_int32_t msw;
    u_int32_t lsw;
  } parts;
  struct
  {
    u_int64_t w;
  } xparts;
} ieee_double_shape_type;

#endif

#if __FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__

typedef union
{
  double value;
  struct
  {
    u_int32_t lsw;
    u_int32_t msw;
  } parts;
  struct
  {
    u_int64_t w;
  } xparts;
} ieee_double_shape_type;

#endif

/* Get two 32 bit ints from a double.  */

#define EXTRACT_WORDS(ix0,ix1,d)				\
do {								\
  ieee_double_shape_type ew_u;					\
  ew_u.value = (d);						\
  (ix0) = ew_u.parts.msw;					\
  (ix1) = ew_u.parts.lsw;					\
} while (0)

/* Get a 64-bit int from a double. */
#define EXTRACT_WORD64(ix,d)					\
do {								\
  ieee_double_shape_type ew_u;					\
  ew_u.value = (d);						\
  (ix) = ew_u.xparts.w;						\
} while (0)

/* Get the more significant 32 bit int from a double.  */

#define GET_HIGH_WORD(i,d)					\
do {								\
  ieee_double_shape_type gh_u;					\
  gh_u.value = (d);						\
  (i) = gh_u.parts.msw;						\
} while (0)

/* Get the less significant 32 bit int from a double.  */

#define GET_LOW_WORD(i,d)					\
do {								\
  ieee_double_shape_type gl_u;					\
  gl_u.value = (d);						\
  (i) = gl_u.parts.lsw;						\
} while (0)

/* Set a double from two 32 bit ints.  */

#define INSERT_WORDS(d,ix0,ix1)					\
do {								\
  ieee_double_shape_type iw_u;					\
  iw_u.parts.msw = (ix0);					\
  iw_u.parts.lsw = (ix1);					\
  (d) = iw_u.value;						\
} while (0)

/* Set a double from a 64-bit int. */
#define INSERT_WORD64(d,ix)					\
do {								\
  ieee_double_shape_type iw_u;					\
  iw_u.xparts.w = (ix);						\
  (d) = iw_u.value;						\
} while (0)

/* Set the more significant 32 bits of a double from an int.  */

#define SET_HIGH_WORD(d,v)					\
do {								\
  ieee_double_shape_type sh_u;					\
  sh_u.value = (d);						\
  sh_u.parts.msw = (v);						\
  (d) = sh_u.value;						\
} while (0)

/* Set the less significant 32 bits of a double from an int.  */

#define SET_LOW_WORD(d,v)					\
do {								\
  ieee_double_shape_type sl_u;					\
  sl_u.value = (d);						\
  sl_u.parts.lsw = (v);						\
  (d) = sl_u.value;						\
} while (0)

/*
 * A union which permits us to convert between a float and a 32 bit
 * int.
 */

typedef union
{
  float value;
  /* FIXME: Assumes 32 bit int.  */
  unsigned int word;
} ieee_float_shape_type;

/* Get a 32 bit int from a float.  */

#define GET_FLOAT_WORD(i,d)					\
do {								\
  ieee_float_shape_type gf_u;					\
  gf_u.value = (d);						\
  (i) = gf_u.word;						\
} while (0)

/* Set a float from a 32 bit int.  */

#define SET_FLOAT_WORD(d,i)					\
do {								\
  ieee_float_shape_type sf_u;					\
  sf_u.word = (i);						\
  (d) = sf_u.value;						\
} while (0)

/* Get expsign as a 16 bit int from a long double.  */

#define	GET_LDBL_EXPSIGN(i,d)					\
do {								\
  union IEEEl2bits ge_u;					\
  ge_u.e = (d);							\
  (i) = ge_u.xbits.expsign;					\
} while (0)

/* Set expsign of a long double from a 16 bit int.  */

#define	SET_LDBL_EXPSIGN(d,v)					\
do {								\
  union IEEEl2bits se_u;					\
  se_u.e = (d);							\
  se_u.xbits.expsign = (v);					\
  (d) = se_u.e;							\
} while (0)


#ifndef __FreeBSD__
#define	STRICT_ASSIGN(type, lval, rval)	((lval) = (rval))
#else
#ifdef FLT_EVAL_METHOD
// Attempt to get strict C99 semantics for assignment with non-C99 compilers.
#if FLT_EVAL_METHOD == 0 || __GNUC__ == 0
#define	STRICT_ASSIGN(type, lval, rval)	((lval) = (rval))
#else
#define	STRICT_ASSIGN(type, lval, rval) do {	\
	volatile type __lval;			\
						\
	if (sizeof(type) >= sizeof(long double))	\
		(lval) = (rval);		\
	else {					\
		__lval = (rval);		\
		(lval) = __lval;		\
	}					\
} while (0)
#endif
#endif
#endif

/*
 * 2sum gives the same result as 2sumF without requiring |a| >= |b| or
 * a == 0, but is slower.
 */
#define	_2sum(a, b) do {	\
	__typeof(a) __s, __w;	\
				\
	__w = (a) + (b);	\
	__s = __w - (a);	\
	(b) = ((a) - (__w - __s)) + ((b) - __s); \
	(a) = __w;		\
} while (0)

/*
 * 2sumF algorithm.
 *
 * "Normalize" the terms in the infinite-precision expression a + b for
 * the sum of 2 floating point values so that b is as small as possible
 * relative to 'a'.  (The resulting 'a' is the value of the expression in
 * the same precision as 'a' and the resulting b is the rounding error.)
 * |a| must be >= |b| or 0, b's type must be no larger than 'a's type, and
 * exponent overflow or underflow must not occur.  This uses a Theorem of
 * Dekker (1971).  See Knuth (1981) 4.2.2 Theorem C.  The name "TwoSum"
 * is apparently due to Skewchuk (1997).
 *
 * For this to always work, assignment of a + b to 'a' must not retain any
 * extra precision in a + b.  This is required by C standards but broken
 * in many compilers.  The brokenness cannot be worked around using
 * STRICT_ASSIGN() like we do elsewhere, since the efficiency of this
 * algorithm would be destroyed by non-null strict assignments.  (The
 * compilers are correct to be broken -- the efficiency of all floating
 * point code calculations would be destroyed similarly if they forced the
 * conversions.)
 *
 * Fortunately, a case that works well can usually be arranged by building
 * any extra precision into the type of 'a' -- 'a' should have type float_t,
 * double_t or long double.  b's type should be no larger than 'a's type.
 * Callers should use these types with scopes as large as possible, to
 * reduce their own extra-precision and efficiciency problems.  In
 * particular, they shouldn't convert back and forth just to call here.
 */
#ifdef DEBUG
#define	_2sumF(a, b) do {				\
	__typeof(a) __w;				\
	volatile __typeof(a) __ia, __ib, __r, __vw;	\
							\
	__ia = (a);					\
	__ib = (b);					\
	assert(__ia == 0 || fabsl(__ia) >= fabsl(__ib));	\
							\
	__w = (a) + (b);				\
	(b) = ((a) - __w) + (b);			\
	(a) = __w;					\
							\
	/* The next 2 assertions are weak if (a) is already long double. */ \
	assert((long double)__ia + __ib == (long double)(a) + (b));	\
	__vw = __ia + __ib;				\
	__r = __ia - __vw;				\
	__r += __ib;					\
	assert(__vw == (a) && __r == (b));		\
} while (0)
#else /* !DEBUG */
#define	_2sumF(a, b) do {	\
	__typeof(a) __w;	\
				\
	__w = (a) + (b);	\
	(b) = ((a) - __w) + (b); \
	(a) = __w;		\
} while (0)
#endif /* DEBUG */

/*
 * Set x += c, where x is represented in extra precision as a + b.
 * x must be sufficiently normalized and sufficiently larger than c,
 * and the result is then sufficiently normalized.
 *
 * The details of ordering are that |a| must be >= |c| (so that (a, c)
 * can be normalized without extra work to swap 'a' with c).  The details of
 * the normalization are that b must be small relative to the normalized 'a'.
 * Normalization of (a, c) makes the normalized c tiny relative to the
 * normalized a, so b remains small relative to 'a' in the result.  However,
 * b need not ever be tiny relative to 'a'.  For example, b might be about
 * 2**20 times smaller than 'a' to give about 20 extra bits of precision.
 * That is usually enough, and adding c (which by normalization is about
 * 2**53 times smaller than a) cannot change b significantly.  However,
 * cancellation of 'a' with c in normalization of (a, c) may reduce 'a'
 * significantly relative to b.  The caller must ensure that significant
 * cancellation doesn't occur, either by having c of the same sign as 'a',
 * or by having |c| a few percent smaller than |a|.  Pre-normalization of
 * (a, b) may help.
 *
 * This is is a variant of an algorithm of Kahan (see Knuth (1981) 4.2.2
 * exercise 19).  We gain considerable efficiency by requiring the terms to
 * be sufficiently normalized and sufficiently increasing.
 */
#define	_3sumF(a, b, c) do {	\
	__typeof(a) __tmp;	\
				\
	__tmp = (c);		\
	_2sumF(__tmp, (a));	\
	(b) += (a);		\
	(a) = __tmp;		\
} while (0)

/*
 * Common routine to process the arguments to nan(), nanf(), and nanl().
 */
void __scan_nan(u_int32_t *__words, int __num_words, const char *__s);

/*
 * Mix 1 or 2 NaNs.  First add 0 to each arg.  This normally just turns
 * signaling NaNs into quiet NaNs by setting a quiet bit.  We do this
 * because we want to never return a signaling NaN, and also because we
 * don't want the quiet bit to affect the result.  Then mix the converted
 * args using addition.  The result is typically the arg whose mantissa
 * bits (considered as in integer) are largest.
 *
 * Technical complications: the result in bits might depend on the precision
 * and/or on compiler optimizations, especially when different register sets
 * are used for different precisions.  Try to make the result not depend on
 * at least the precision by always doing the main mixing step in long double
 * precision.  Try to reduce dependencies on optimizations by adding the
 * the 0's in different precisions (unless everything is in long double
 * precision).
 */
#define	nan_mix(x, y)	(((x) + 0.0L) + ((y) + 0))

#ifdef __GNUCLIKE_ASM

/* Asm versions of some functions. */

#ifdef __amd64__
static __inline int
irint(double x)
{
	int n;

	__asm__("cvtsd2si %1,%0" : "=r" (n) : "x" (x));
	return (n);
}
#define	HAVE_EFFICIENT_IRINT
#endif

#ifdef __i386__
static __inline int
irint(double x)
{
	int n;

	__asm__("fistl %0" : "=m" (n) : "t" (x));
	return (n);
}
#define	HAVE_EFFICIENT_IRINT
#endif

#endif /* __GNUCLIKE_ASM */

/*
 * ieee style elementary functions
 *
 * We rename functions here to improve other sources' diffability
 * against fdlibm.
 */
#define	__ieee754_sqrt	sqrt
#define	__ieee754_acos	acos
#define	__ieee754_acosh	acosh
#define	__ieee754_log	log
#define	__ieee754_log2	log2
#define	__ieee754_atanh	atanh
#define	__ieee754_asin	asin
#define	__ieee754_atan2	atan2
#define	__ieee754_exp	exp
#define	__ieee754_cosh	cosh
#define	__ieee754_fmod	fmod
#define	__ieee754_pow	pow
#define	__ieee754_lgamma lgamma
#define	__ieee754_lgamma_r lgamma_r
#define	__ieee754_log10	log10
#define	__ieee754_sinh	sinh
#define	__ieee754_hypot	hypot
#define	__ieee754_j0	j0
#define	__ieee754_j1	j1
#define	__ieee754_y0	y0
#define	__ieee754_y1	y1
#define	__ieee754_jn	jn
#define	__ieee754_yn	yn
#define	__ieee754_remainder remainder
#define	__ieee754_sqrtf	sqrtf
#define	__ieee754_acosf	acosf
#define	__ieee754_acoshf acoshf
#define	__ieee754_logf	logf
#define	__ieee754_atanhf atanhf
#define	__ieee754_asinf	asinf
#define	__ieee754_atan2f atan2f
#define	__ieee754_expf	expf
#define	__ieee754_coshf	coshf
#define	__ieee754_fmodf	fmodf
#define	__ieee754_powf	powf
#define	__ieee754_lgammaf lgammaf
#define	__ieee754_lgammaf_r lgammaf_r
#define	__ieee754_log10f log10f
#define	__ieee754_log2f log2f
#define	__ieee754_sinhf	sinhf
#define	__ieee754_hypotf hypotf
#define	__ieee754_j0f	j0f
#define	__ieee754_j1f	j1f
#define	__ieee754_y0f	y0f
#define	__ieee754_y1f	y1f
#define	__ieee754_jnf	jnf
#define	__ieee754_ynf	ynf
#define	__ieee754_remainderf remainderf

/* fdlibm kernel function */
int	__kernel_rem_pio2(double*,double*,int,int,int);

/* double precision kernel functions */
#ifdef INLINE_REM_PIO2
__inline
#endif
int	__ieee754_rem_pio2(double,double*);
double	__kernel_sin(double,double,int);
double	__kernel_cos(double,double);
double	__kernel_tan(double,double,int);
double	__ldexp_exp(double,int);
double complex __ldexp_cexp(double complex,int);

/* float precision kernel functions */
#ifndef INLINE_REM_PIO2F
int     __ieee754_rem_pio2f(float,double*);
#endif
#ifndef INLINE_KERNEL_SINDF
float   __kernel_sindf(double);
#endif
#ifndef INLINE_KERNEL_COSDF
float   __kernel_cosdf(double);
#endif
#ifndef INLINE_KERNEL_TANDF
float   __kernel_tandf(double,int);
#endif
float   __ldexp_expf(float,int);
#ifdef _COMPLEX_H
float complex __ldexp_cexpf(float complex,int);
#endif

/* long double precision kernel functions */
long double __kernel_sinl(long double, long double, int);
long double __kernel_cosl(long double, long double);
long double __kernel_tanl(long double, long double, int);

#endif /* !_MATH_PRIVATE_H_ */
