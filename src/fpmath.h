/*-
 * Copyright (c) 2003 Mike Barcroft <mike@FreeBSD.org>
 * Copyright (c) 2002 David Schultz <das@FreeBSD.ORG>
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
 * $FreeBSD: src/lib/libc/include/fpmath.h,v 1.4 2008/12/23 22:20:59 marcel Exp $
 */
#ifndef _FPMATH_H_
#define _FPMATH_H_

#if defined(__aarch64__)
#include "aarch64_fpmath.h"
#elif defined(__i386__) || defined(__x86_64__)
#ifdef __LP64__
#include "amd64_fpmath.h"
#else 
#include "i386_fpmath.h"
#endif
#elif defined(__powerpc__)
#include "powerpc_fpmath.h"
#elif defined(__mips__)
#include "mips_fpmath.h"
#elif defined(__s390__)
#include "s390_fpmath.h"
#elif defined(__riscv)
#include "riscv_fpmath.h"
#endif

/* Definitions provided directly by GCC and Clang. */
#if !(defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__))

#if defined(__GLIBC__)

#include <features.h>
#include <endian.h>
#define __ORDER_LITTLE_ENDIAN__  __LITTLE_ENDIAN
#define __ORDER_BIG_ENDIAN__     __BIG_ENDIAN
#define __BYTE_ORDER__           __BYTE_ORDER

#elif defined(__APPLE__)

#include <machine/endian.h>
#define __ORDER_LITTLE_ENDIAN__  LITTLE_ENDIAN
#define __ORDER_BIG_ENDIAN__     BIG_ENDIAN
#define __BYTE_ORDER__           BYTE_ORDER

#elif defined(_WIN32)

#define __ORDER_LITTLE_ENDIAN__  1234
#define __ORDER_BIG_ENDIAN__     4321
#define __BYTE_ORDER__           __ORDER_LITTLE_ENDIAN__

#endif

#endif /* __BYTE_ORDER__, __ORDER_LITTLE_ENDIAN__ and __ORDER_BIG_ENDIAN__ */

#ifndef __FLOAT_WORD_ORDER__
#define __FLOAT_WORD_ORDER__     __BYTE_ORDER__
#endif

union IEEEf2bits {
	float	f;
	struct {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
		unsigned int	man	:23;
		unsigned int	exp	:8;
		unsigned int	sign	:1;
#else /* _BIG_ENDIAN */
		unsigned int	sign	:1;
		unsigned int	exp	:8;
		unsigned int	man	:23;
#endif
	} bits;
};

#define	DBL_MANH_SIZE	20
#define	DBL_MANL_SIZE	32

union IEEEd2bits {
	double	d;
	struct {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#if __FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__
		unsigned int	manl	:32;
#endif
		unsigned int	manh	:20;
		unsigned int	exp	:11;
		unsigned int	sign	:1;
#if __FLOAT_WORD_ORDER__ == __ORDER_BIG_ENDIAN__
		unsigned int	manl	:32;
#endif
#else /* _BIG_ENDIAN */
		unsigned int	sign	:1;
		unsigned int	exp	:11;
		unsigned int	manh	:20;
		unsigned int	manl	:32;
#endif
	} bits;
};

#endif
