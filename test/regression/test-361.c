/*
 * Regression test for issue #361: nextafter/nexttoward must return signed
 * infinity at the finite endpoints in every rounding mode, while raising
 * overflow and inexact (C Annex F.10.8.3-4).
 */
#include <float.h>
#include <openlibm_fenv.h>
#include <openlibm_math.h>

#include "regress-util.h"

#define CHECK_ENDPOINT(expr) do { \
	long double result; \
	int flags; \
	CHECK(feclearexcept(FE_ALL_EXCEPT) == 0); \
	result = (expr); \
	flags = fetestexcept(FE_ALL_EXCEPT); \
	CHECK(isinf(result)); \
	CHECK(!!signbit(result) == (sign < 0)); \
	CHECK((flags & (FE_OVERFLOW | FE_INEXACT)) == \
	    (FE_OVERFLOW | FE_INEXACT)); \
	CHECK(fegetround() == modes[m]); \
} while (0)

int
main(void)
{
	const int modes[] = {
		FE_TONEAREST, FE_UPWARD, FE_DOWNWARD, FE_TOWARDZERO
	};
	int saved_round = fegetround();
	unsigned int m;
	int sign;

	for (m = 0; m < sizeof(modes) / sizeof(modes[0]); m++) {
		CHECK(fesetround(modes[m]) == 0);
		for (sign = -1; sign <= 1; sign += 2) {
			CHECK_ENDPOINT(nextafter(sign * DBL_MAX,
			    sign * INFINITY));
			CHECK_ENDPOINT(nextafterf(sign * FLT_MAX,
			    sign * INFINITY));
#ifdef OLM_LONG_DOUBLE
			CHECK_ENDPOINT(nexttowardf(sign * FLT_MAX,
			    sign * (long double)INFINITY));
/* Extended long double routines are built only on x86 and aarch64. */
#if LDBL_MANT_DIG == 53 || LDBL_MANT_DIG == 64 || defined(__aarch64__)
			CHECK_ENDPOINT(nextafterl(sign * LDBL_MAX,
			    sign * (long double)INFINITY));
			CHECK_ENDPOINT(nexttoward(sign * DBL_MAX,
			    sign * (long double)INFINITY));
			CHECK_ENDPOINT(nexttowardl(sign * LDBL_MAX,
			    sign * (long double)INFINITY));
#endif
#endif
		}
	}
	CHECK(fesetround(saved_round) == 0);
	return 0;
}
