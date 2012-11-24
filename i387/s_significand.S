/*
 * Written by:
 *	J.T. Conklin (jtc@netbsd.org)
 *	Public domain.
 */

#include <i387/bsd_asm.h>
//__FBSDID("$FreeBSD: src/lib/msun/i387/s_significand.S,v 1.10 2011/01/07 16:13:12 kib Exp $")

ENTRY(significand)
	fldl	4(%esp)
	fxtract
	fstp	%st(1)
	ret
END(significand)

	
