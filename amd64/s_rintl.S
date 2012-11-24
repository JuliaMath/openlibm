/*
 * Written by:
 *	J.T. Conklin (jtc@netbsd.org)
 *	Public domain.
 */

#include <amd64/bsd_asm.h>

ENTRY(rintl)
	fldt	8(%rsp)
	frndint
	ret

	
