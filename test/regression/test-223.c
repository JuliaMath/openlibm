/*
 * Regression test for issue #223: tgammal spuriously overflows for
 * large-ish negative non-integer arguments.
 * https://github.com/JuliaMath/openlibm/issues/223
 *
 * For a large negative non-integer x, |tgamma(x)| is tiny (the reflection
 * formula divides pi by the enormous tgamma(1-x)).  The ld80 implementation
 * applied its positive-overflow threshold (MAXGAML) to negative arguments
 * too and, just below that threshold, let z * tgamma(q) overflow to Inf
 * before the final division -- both paths wrongly produced +Inf (or a
 * flush-to-zero) instead of the correct tiny subnormal.
 *
 * glibc on x86-64 (the correctness oracle) returns the correctly-rounded
 * subnormal 0x0.01dbd551da54538p-16385L for the argument below; pin that.
 */
#include <float.h>

#include "regress-util.h"

/*
 * ld80-specific (the fix is in ld80/e_tgammal.c, and the exact hex below holds
 * only for the 80-bit format).  openlibm provides long double routines only on
 * x86 and (linux) aarch64, so every other format must compile to a bare skip
 * with no reference to tgammal (it may be absent, which would break the link).
 */
#if LDBL_MANT_DIG != 64

int
main(void)
{
	return REGRESS_SKIP;
}

#else

#include <math.h>

int
main(void)
{
	long double x = -0xd.b6e8f5c28f5c29p+7L;	/* approx -1755.455 */
	long double y = tgammal(x);

	/* Must be a finite, positive, tiny value -- not Inf and not zero. */
	CHECK(isfinite(y));
	CHECK(y > 0.0L);
	CHECK(y < 1.0L);

	/* Correctly-rounded result observed from glibc tgammal on x86-64. */
	CHECK(y == 0x0.01dbd551da54538p-16385L);

	/* The positive-overflow path (moved out of the shared code) must still
	   return +Inf for large positive x. */
	CHECK(isinf(tgammal(2000.0L)) == 1);

	return 0;
}

#endif
