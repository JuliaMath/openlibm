/*
 * Regression test for issue #306: expm1l overflows internally and returns NaN
 * for large arguments that are still below the overflow threshold.
 * https://github.com/JuliaMath/openlibm/issues/306
 *
 * For x just below MAXLOG the true result is finite (near LDBL_MAX), but the
 * old reconstruction formed 2^k (with k == LDBL_MAX_EXP) which overflowed to
 * +Inf; the subsequent (+Inf) + (-Inf) yielded NaN.  The fix scales the
 * exponent of exp(remainder) directly in that regime.
 *
 * The expected hex values below are bit-for-bit what the system glibc expm1l
 * returns on x86-64 (LDBL_MANT_DIG == 64).  Only meaningful for the 80-bit
 * long double format, so skip elsewhere.
 */
#include <float.h>
#include <math.h>

#include "regress-util.h"

int
main(void)
{
	/*
	 * Target the x86 80-bit long double path only: that is the code this fix
	 * changes, and the only format for which the exact hex value below holds.
	 * openlibm provides long double routines only on x86 and (linux) aarch64,
	 * so every other format must compile to a bare skip with no reference to
	 * expm1l (it may be absent from the library, breaking the link).
	 */
#if LDBL_MANT_DIG != 64
	return REGRESS_SKIP;
#else
	long double x = 0x2.c5c85fdf170c604cp+12l;	/* ~11346.6, just below threshold */
	long double y = expm1l(x);

	/* Must be a large finite value, not NaN/Inf.  isnan() is type-generic;
	   isnanl is a glibc-only symbol that would not link on macOS. */
	CHECK(!isnan(y));
	CHECK(isfinite(y));
	CHECK(y > 0.0L);

	/* Exact value matching glibc. */
	CHECK(y == 0xf.ffffcfce79e56d5p+16380l);

	/* A slightly larger argument must still overflow to +Inf. */
	CHECK(isinf(expm1l(11357.0L)) == 1);

	return 0;
#endif
}
