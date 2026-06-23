/*-
 * Copyright (c) 2005-2011 David Schultz <das@FreeBSD.ORG>
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
 */

#include "cdefs-compat.h"
//__FBSDID("$FreeBSD: src/lib/msun/src/s_fmaf.c,v 1.3 2011/10/15 04:16:58 das Exp $");

#include <openlibm_fenv.h>
#include <openlibm_math.h>

#include "math_private.h"

/*
 * Fused multiply-add: Compute x * y + z with a single rounding error.
 *
 * A double has more than twice as much precision than a float, so
 * direct double-precision arithmetic suffices, except where double
 * rounding occurs.
 */
OLM_DLLEXPORT float
fmaf(float x, float y, float z)
{
	double xy, result;
	u_int32_t hr, lr;

	xy = (double)x * y;
	result = xy + z;
	EXTRACT_WORDS(hr, lr, result);

	/*
	 * Compute the exact rounding error of the sum xy + z using Knuth's
	 * 2Sum.  result is the exact value of xy + z if and only if err == 0.
	 *
	 * The naive shortcut "result - xy == z" is *not* a valid exactness
	 * test: when |z| > |xy| the subtraction result - xy is itself rounded
	 * and can return z even though result lost low-order bits of the true
	 * sum (e.g. fmaf(0.9474001f, 4.639901e-7f, -0.24325085f)).  Trusting
	 * it skips the correction below and lets the inexact double result be
	 * rounded a second time to float, producing a 1-ulp double-rounding
	 * error.  2Sum is order-independent, so it avoids that trap.
	 *
	 * The volatile barriers force each step to round in double precision
	 * so the decomposition stays exact even under FP contraction.
	 */
	volatile double zz = result - xy;
	volatile double rz = result - zz;
	volatile double tz = z - zz;
	double err = (xy - rz) + tz;

	/* Common case: The double precision result is fine. */
	if ((lr & 0x1fffffff) != 0x10000000 ||	/* not a halfway case */
	    (hr & 0x7ff00000) == 0x7ff00000 ||	/* NaN */
	    err == 0 ||				/* exact */
	    fegetround() != FE_TONEAREST)	/* not round-to-nearest */
		return (result);

	/*
	 * If result is inexact, and exactly halfway between two float values,
	 * we need to adjust the low-order bit in the direction of the error.
	 */
	fesetround(FE_TOWARDZERO);
	volatile double vxy = xy;  /* XXX work around gcc CSE bug */
	double adjusted_result = vxy + z;
	fesetround(FE_TONEAREST);
	if (result == adjusted_result)
		SET_LOW_WORD(adjusted_result, lr + 1);
	return (adjusted_result);
}
