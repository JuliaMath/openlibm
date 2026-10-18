/*
 * scalblnl() must clamp large long exponents before converting to int.
 * On LP64, truncating +/-2^32 to int used to scale by zero instead.
 */
#include <float.h>
#include <limits.h>
#include <openlibm_fenv.h>
#include <openlibm_math.h>

#include "regress-util.h"

int
main(void)
{
#if LONG_MAX <= INT_MAX || LDBL_MANT_DIG != 64
	return REGRESS_SKIP;
#else
	const long exponents[] = {
		4294967296L, -4294967296L, LONG_MAX, LONG_MIN
	};
	long double actual;
	unsigned i;
	int flags, expected;

	CHECK(fesetround(FE_TONEAREST) == 0);
	CHECK(scalblnl(1.0L, 4) == 16.0L);
	CHECK(scalblnl(1.0L, -4) == 0.0625L);
	for (i = 0; i < sizeof(exponents) / sizeof(exponents[0]); i++) {
		CHECK(feclearexcept(FE_ALL_EXCEPT) == 0);
		actual = scalblnl(1.0L, exponents[i]);
		flags = fetestexcept(FE_ALL_EXCEPT);
		if (exponents[i] > 0) {
			CHECK(isinf(actual) && actual > 0.0L);
			expected = FE_OVERFLOW | FE_INEXACT;
		} else {
			CHECK(actual == 0.0L && !signbit(actual));
			expected = FE_UNDERFLOW | FE_INEXACT;
		}
		CHECK((flags & expected) == expected);
	}
	return 0;
#endif
}
