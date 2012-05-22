/*-
 * Copyright (c) 2004 David Schultz <das@FreeBSD.ORG>
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
 */

#include <openlibm.h>

#include "fpmath.h"

//#define FP_INFINITE     0x01
//#define FP_NAN          0x02
//#define FP_NORMAL       0x04
//#define FP_SUBNORMAL    0x08
//#define FP_ZERO         0x10
//#define fpclassify(x) \
//    ((sizeof (x) == sizeof (float)) ? __fpclassifyf(x) \
//    : (sizeof (x) == sizeof (double)) ? __fpclassifyd(x) \
//    : __fpclassifyl(x))
//

int
__fpclassifyd(double d)
{
	union IEEEd2bits u;

	u.d = d;
	if (u.bits.exp == 2047) {
		if (u.bits.manl == 0 && u.bits.manh == 0) {
			return FP_INFINITE;
		} else {
			return FP_NAN;
		}
	} else if (u.bits.exp != 0) {
		return FP_NORMAL;
	} else if (u.bits.manl == 0 && u.bits.manh == 0) {
		return FP_ZERO;
	} else {
		return FP_SUBNORMAL;
	}
}
		

int
__fpclassifyf(float f)
{
	union IEEEf2bits u;

	u.f = f;
	if (u.bits.exp == 255) {
	   if (u.bits.man == 0) {
		   return FP_INFINITE;
	   } else {
		   return FP_NAN;
	   }
	} else if (u.bits.exp != 0) {
		return FP_NORMAL;
	} else if (u.bits.man == 0) {
		return FP_ZERO;
	} else {
		return FP_SUBNORMAL;
	}
}

int
__fpclassifyl(long double e)
{
	union IEEEl2bits u;

	u.e = e;
	mask_nbit_l(u);
	if (u.bits.exp == 32767) {
		if (u.bits.manl == 0 && u.bits.manh == 0) {
			return FP_INFINITE;
		} else {
			return FP_NAN;
		}
	} else if (u.bits.exp != 0) {
		return FP_NORMAL;
	} else if (u.bits.manl == 0 && u.bits.manh == 0) {
		return FP_ZERO;
	} else {
		return FP_SUBNORMAL;
	}
}


