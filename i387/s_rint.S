/*
 * Written by:
 *	J.T. Conklin (jtc@netbsd.org)
 *	Public domain.
 */

#include <i387/bsd_asm.h>
//__FBSDID("$FreeBSD: src/lib/msun/i387/s_rint.S,v 1.9 2011/01/07 16:13:12 kib Exp $")

ENTRY(rint)
	fldl	4(%esp)
	frndint
	ret
END(rint)

	
