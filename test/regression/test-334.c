/*
 * Regression test for issue #334: powl() returns NaN instead of +0 on
 * extreme underflow.
 * https://github.com/JuliaMath/openlibm/issues/334
 *
 * For very large |y * log2(x)|, the extended-precision reduction step in
 * ld80/e_powl.c (reducl(G), which computes ldexpl(G, LNXT)) overflowed to
 * +-Inf.  The following Inf-Inf / Inf+(-Inf) arithmetic produced NaN, which
 * defeated the w > MEXP / w < MNEXP over/underflow tests and made powl()
 * return NaN.  The mathematically correct result here is +0 (underflow),
 * which is what glibc returns.
 *
 * This manifests only when long double has more precision than double (the
 * x86 80-bit / ld80 path); otherwise powl == pow and the bug cannot occur.
 */
#include <float.h>
#include <math.h>

#include "regress-util.h"

int
main(void)
{
	/*
	 * ld80-specific (the fix is in ld80/e_powl.c).  openlibm provides long
	 * double routines only on x86 and (linux) aarch64, so restrict to the
	 * 80-bit format; other arches compile to a bare skip with no powl
	 * reference (powl may be absent, which would break the link).
	 */
#if LDBL_MANT_DIG != 64
	return REGRESS_SKIP;
#else
	long double x = 0x1.98p-3072l;
	long double y = 0xa.ab43b6dba9a6383p+16364l;
	long double z = powl(x, y);

	/* Must underflow to +0, not return NaN. */
	CHECK(!isnan(z));
	CHECK(z == 0.0L);
	CHECK(!signbit(z));   /* result is +0, not -0 */

	/* Exercise the symmetric early-overflow guard: a result far beyond the
	   range (here 10^1e6) must return +Inf, matching glibc. */
	CHECK(isinf(powl(10.0L, 1000000.0L)) == 1);
	return 0;
#endif
}
