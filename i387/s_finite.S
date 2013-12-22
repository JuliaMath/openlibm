/*
 * Written by:
 *	J.T. Conklin (jtc@netbsd.org)
 *	Public domain.
 */

#include <i387/bsd_asm.h>
//__FBSDID("$FreeBSD: src/lib/msun/i387/s_finite.S,v 1.10 2011/01/07 16:13:12 kib Exp $")

ENTRY(finite)
	movl	8(%esp),%eax
	andl	$0x7ff00000, %eax
	cmpl	$0x7ff00000, %eax
	setneb	%al
	andl	$0x000000ff, %eax
	ret
END(finite)

	
/* Enable stack protection */
#if defined(__linux__) && defined(__ELF__)
.section .note.GNU-stack,"",%progbits
#endif
