/*
 * Regression test for issue #285: ld128 coshl returned 1+x instead of 1 for
 * tiny |x|.
 * https://github.com/JuliaMath/openlibm/issues/285
 *
 * The early-exit branch in ld128/e_coshl.c returned `1 + expm1l(x)` (~= 1+x)
 * for tiny representable x instead of exactly 1.0. cosh(x) = 1 + x^2/2 + O(x^4);
 * for these inputs the x^2/2 term is below half an ulp at 1.0 (2^-113 in
 * 113-bit precision), so the correctly rounded result is exactly 1.0.
 *
 * This only manifests where long double is 128-bit (LDBL_MANT_DIG == 113),
 * e.g. aarch64/riscv64 Linux. Elsewhere (80-bit x87, 64-bit double) the test
 * is skipped.
 */
#include <float.h>

#include "regress-util.h"

/*
 * The fix is in ld128/e_coshl.c, which openlibm compiles only on (linux)
 * aarch64.  Restrict the test to that target at compile time so that on every
 * other arch it becomes a bare skip with no reference to coshl (which may be
 * absent from the library and would otherwise break the link).
 */
#if defined(__aarch64__) && LDBL_MANT_DIG == 113

#include <math.h>

int
main(void)
{
	int k;

	/*
	 * For k in [57,128], coshl(2^-k) rounds to exactly 1.0L in 113-bit
	 * precision. The pre-fix code returned 1+tiny for k > 71.
	 */
	for (k = 57; k <= 128; k++)
		CHECK(coshl(ldexpl(1.0L, -k)) == 1.0L);

	return 0;
}

#else

int
main(void)
{
	return REGRESS_SKIP;
}

#endif
