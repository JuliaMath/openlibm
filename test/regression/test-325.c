/*
 * Regression test for issue #325: tgamma() error signaling.
 * https://github.com/JuliaMath/openlibm/issues/325
 *
 * openlibm reports math errors via FP exception flags, not errno
 * (math_errhandling == MATH_ERREXCEPT).  tgamma already raised FE_DIVBYZERO at
 * the pole (x == +-0) and FE_INVALID for negative-integer (domain) arguments,
 * but overflow for large x went through `x / zero`, which returns +Inf while
 * wrongly raising FE_DIVBYZERO instead of FE_OVERFLOW.  This pins the correct,
 * glibc-matching exception for every error class.
 *
 * FP exception flags are emulated unreliably by qemu-user on some arches (see
 * #347), so the test first checks that the environment tracks a plain
 * divide-by-zero and skips if it does not.
 */
#include <math.h>
#include <openlibm_fenv.h>

#include "regress-util.h"

static int
fenv_tracks_flags(void)
{
	volatile double z = 0.0, r;

	feclearexcept(FE_ALL_EXCEPT);
	r = 1.0 / z;			/* a known divide-by-zero */
	(void)r;
	return fetestexcept(FE_DIVBYZERO) != 0;
}

static int
raises(double x, int excepts)
{
	feclearexcept(FE_ALL_EXCEPT);
	(void)tgamma(x);
	return fetestexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW) == excepts;
}

int
main(void)
{
	if (!fenv_tracks_flags())
		return REGRESS_SKIP;

	/* Pole at x == +-0: +-Inf, FE_DIVBYZERO. */
	CHECK(isinf(tgamma(0.0)) == 1);
	CHECK(raises(0.0, FE_DIVBYZERO));
	CHECK(isinf(tgamma(-0.0)) == -1);
	CHECK(raises(-0.0, FE_DIVBYZERO));

	/* Negative integer (domain error): NaN, FE_INVALID. */
	CHECK(isnan(tgamma(-2.0)));
	CHECK(raises(-2.0, FE_INVALID));

	/* Overflow for large finite x: +Inf, FE_OVERFLOW (was FE_DIVBYZERO). */
	CHECK(isinf(tgamma(200.0)) == 1);
	CHECK(raises(200.0, FE_OVERFLOW));

	/* In-range result raises none of these. */
	CHECK(tgamma(5.0) == 24.0);
	CHECK(raises(5.0, 0));

	return 0;
}
