/*-
 * Copyright (c) 2004-2005 David Schultz <das@FreeBSD.ORG>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * $FreeBSD: src/lib/msun/i387/fenv.h,v 1.8 2011/10/10 15:43:09 das Exp $
 */

#ifndef	_FENV_H_
#define	_FENV_H_

#include "openlibm_defs.h"
#include "cdefs-compat.h"
#include "types-compat.h"

#ifndef	__fenv_static
#define	__fenv_static	static
#endif

/*                   
 * To preserve binary compatibility with FreeBSD 5.3, we pack the
 * mxcsr into some reserved fields, rather than changing sizeof(fenv_t).
 */
typedef struct {
	uint16_t	__control;
	uint16_t      __mxcsr_hi;
	uint16_t	__status;
	uint16_t      __mxcsr_lo;
	uint32_t	__tag;
	char		__other[16];
} fenv_t;

#define	__get_mxcsr(env)	(((env).__mxcsr_hi << 16) |	\
				 ((env).__mxcsr_lo))
#define	__set_mxcsr(env, x)	do {				\
	(env).__mxcsr_hi = (uint32_t)(x) >> 16;		\
	(env).__mxcsr_lo = (uint16_t)(x);			\
} while (0)

typedef	uint16_t	fexcept_t;

/* Exception flags */
#define	FE_INVALID	0x01
#define	FE_DENORMAL	0x02
#define	FE_DIVBYZERO	0x04
#define	FE_OVERFLOW	0x08
#define	FE_UNDERFLOW	0x10
#define	FE_INEXACT	0x20
#define	FE_ALL_EXCEPT	(FE_DIVBYZERO | FE_DENORMAL | FE_INEXACT | \
			 FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW)

/* Rounding modes */
#define	FE_TONEAREST	0x0000
#define	FE_DOWNWARD	0x0400
#define	FE_UPWARD	0x0800
#define	FE_TOWARDZERO	0x0c00
#define	_ROUND_MASK	(FE_TONEAREST | FE_DOWNWARD | \
			 FE_UPWARD | FE_TOWARDZERO)

/*
 * As compared to the x87 control word, the SSE unit's control word
 * has the rounding control bits offset by 3 and the exception mask
 * bits offset by 7.
 */
#define	_SSE_ROUND_SHIFT	3
#define	_SSE_EMASK_SHIFT	7

__BEGIN_DECLS

/* After testing for SSE support once, we cache the result in __has_sse. */
enum __sse_support { __SSE_YES, __SSE_NO, __SSE_UNK };
OLM_DLLEXPORT extern enum __sse_support __has_sse;
OLM_DLLEXPORT int __test_sse(void);
#ifdef __SSE__
#define	__HAS_SSE()	1
#else
#define	__HAS_SSE()	(__has_sse == __SSE_YES ||			\
			 (__has_sse == __SSE_UNK && __test_sse()))
#endif

/* Default floating-point environment */
OLM_DLLEXPORT extern const fenv_t	__fe_dfl_env;
#define	FE_DFL_ENV	(&__fe_dfl_env)

#define	__fldcw(__cw)		__asm __volatile("fldcw %0" : : "m" (__cw))
#define	__fldenv(__env)		__asm __volatile("fldenv %0" : : "m" (__env))
#define	__fldenvx(__env)	__asm __volatile("fldenv %0" : : "m" (__env)  \
				: "st", "st(1)", "st(2)", "st(3)", "st(4)",   \
				"st(5)", "st(6)", "st(7)")
#define	__fnclex()		__asm __volatile("fnclex")
#define	__fnstenv(__env)	__asm __volatile("fnstenv %0" : "=m" (*(__env)))
#define	__fnstcw(__cw)		__asm __volatile("fnstcw %0" : "=m" (*(__cw)))
#define	__fnstsw(__sw)		__asm __volatile("fnstsw %0" : "=am" (*(__sw)))
#define	__fwait()		__asm __volatile("fwait")
#define	__ldmxcsr(__csr)	__asm __volatile("ldmxcsr %0" : : "m" (__csr))
#define	__stmxcsr(__csr)	__asm __volatile("stmxcsr %0" : "=m" (*(__csr)))

__fenv_static inline int
feclearexcept(int __excepts)
{
	fenv_t __env;
	uint32_t __mxcsr;

	if (__excepts == FE_ALL_EXCEPT) {
		__fnclex();
	} else {
		__fnstenv(&__env);
		__env.__status &= ~__excepts;
		__fldenv(__env);
	}
	if (__HAS_SSE()) {
		__stmxcsr(&__mxcsr);
		__mxcsr &= ~__excepts;
		__ldmxcsr(__mxcsr);
	}
	return (0);
}

__fenv_static inline int
fegetexceptflag(fexcept_t *__flagp, int __excepts)
{
	uint32_t __mxcsr;
	uint16_t __status;

	__fnstsw(&__status);
	if (__HAS_SSE())
		__stmxcsr(&__mxcsr);
	else
		__mxcsr = 0;
	*__flagp = (__mxcsr | __status) & __excepts;
	return (0);
}

OLM_DLLEXPORT int fesetexceptflag(const fexcept_t *__flagp, int __excepts);
OLM_DLLEXPORT int feraiseexcept(int __excepts);

__fenv_static inline int
fetestexcept(int __excepts)
{
	uint32_t __mxcsr;
	uint16_t __status;

	__fnstsw(&__status);
	if (__HAS_SSE())
		__stmxcsr(&__mxcsr);
	else
		__mxcsr = 0;
	return ((__status | __mxcsr) & __excepts);
}

__fenv_static inline int
fegetround(void)
{
	uint16_t __control;

	/*
	 * We assume that the x87 and the SSE unit agree on the
	 * rounding mode.  Reading the control word on the x87 turns
	 * out to be about 5 times faster than reading it on the SSE
	 * unit on an Opteron 244.
	 */
	__fnstcw(&__control);
	return (__control & _ROUND_MASK);
}

__fenv_static inline int
fesetround(int __round)
{
	uint32_t __mxcsr;
	uint16_t __control;

	if (__round & ~_ROUND_MASK)
		return (-1);

	__fnstcw(&__control);
	__control &= ~_ROUND_MASK;
	__control |= __round;
	__fldcw(__control);

	if (__HAS_SSE()) {
		__stmxcsr(&__mxcsr);
		__mxcsr &= ~(_ROUND_MASK << _SSE_ROUND_SHIFT);
		__mxcsr |= __round << _SSE_ROUND_SHIFT;
		__ldmxcsr(__mxcsr);
	}

	return (0);
}

OLM_DLLEXPORT int fegetenv(fenv_t *__envp);
OLM_DLLEXPORT int feholdexcept(fenv_t *__envp);

__fenv_static inline int
fesetenv(const fenv_t *__envp)
{
	fenv_t __env = *__envp;
	uint32_t __mxcsr;

	__mxcsr = __get_mxcsr(__env);
	__set_mxcsr(__env, 0xffffffff);
	/*
	 * XXX Using fldenvx() instead of fldenv() tells the compiler that this
	 * instruction clobbers the i387 register stack.  This happens because
	 * we restore the tag word from the saved environment.  Normally, this
	 * would happen anyway and we wouldn't care, because the ABI allows
	 * function calls to clobber the i387 regs.  However, fesetenv() is
	 * inlined, so we need to be more careful.
	 */
	__fldenvx(__env);
	if (__HAS_SSE())
		__ldmxcsr(__mxcsr);
	return (0);
}

OLM_DLLEXPORT int feupdateenv(const fenv_t *__envp);

#if __BSD_VISIBLE

OLM_DLLEXPORT int feenableexcept(int __mask);
OLM_DLLEXPORT int fedisableexcept(int __mask);

/* We currently provide no external definition of fegetexcept(). */
static inline int
fegetexcept(void)
{
	uint16_t __control;

	/*
	 * We assume that the masks for the x87 and the SSE unit are
	 * the same.
	 */
	__fnstcw(&__control);
	return (~__control & FE_ALL_EXCEPT);
}

#endif /* __BSD_VISIBLE */

__END_DECLS

#endif	/* !_FENV_H_ */
