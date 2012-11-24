/*
 * Written by:
 *	J.T. Conklin (jtc@netbsd.org)
 *	Public domain.
 */

#include <amd64/bsd_asm.h>
//__FBSDID("$FreeBSD: src/lib/msun/amd64/s_logbl.S,v 1.4 2011/01/07 16:13:12 kib Exp $")

ENTRY(logbl)
	fldt	8(%rsp)
	fxtract
	fstp	%st
	ret

	
