/*
 * Based on the i387 version written by:
 *	J.T. Conklin (jtc@netbsd.org)
 *	Public domain.
 */

#include <amd64/bsd_asm.h>

//__FBSDID("$FreeBSD: src/lib/msun/amd64/e_remainderl.S,v 1.2 2011/01/07 16:13:12 kib Exp $")

ENTRY(remainderl)
	fldt	24(%rsp)
	fldt	8(%rsp)
1:	fprem1
	fstsw	%ax
	testw	$0x400,%ax
	jne	1b
	fstp	%st(1)
	ret

	
