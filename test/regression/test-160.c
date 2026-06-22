/*
 * Regression test for issue #160: fmaf double-rounding error.
 * https://github.com/JuliaMath/openlibm/issues/160
 *
 * fmaf must compute x*y + z with a *single* rounding to float.  The old
 * implementation evaluated (double)x*y + z in double precision and rounded
 * to float, using "result - xy == z" as an exactness shortcut to decide
 * whether a halfway result needed correction.  That test is bogus when
 * |z| > |xy|: the subtraction result - xy is itself rounded and can return
 * z even though result already lost low-order bits of the true sum.  The
 * shortcut then skipped the correction and let the inexact double result be
 * rounded a second time to float, giving a 1-ulp double-rounding error.
 *
 * For fmaf(0.9474001f, 4.639901e-7f, -0.24325085f) the correctly single-
 * rounded result is -0x1.f22d46p-3f (== -0.24325041f); the buggy code
 * returned -0x1.f22d44p-3f (== -0.2432504f).
 *
 * The fix replaces the shortcut with an exact Knuth 2Sum error term.  The
 * spot-check values below were verified against a hardware (single-rounded)
 * FMA.
 */
#include <math.h>
#include <openlibm_fenv.h>	/* openlibm's fegetround (inline on x86) */

#include "regress-util.h"

/*
 * Call fmaf through volatile operands so the compiler cannot fold the call
 * at compile time (which would evaluate it with its own arithmetic and hide
 * a bug in the library).  This exercises the real openlibm fmaf.
 */
static float
fma3(float a, float b, float c)
{
	volatile float va = a, vb = b, vc = c;
	return fmaf(va, vb, vc);
}

int
main(void)
{
	/*
	 * fmaf only applies its halfway-case correction in round-to-nearest,
	 * which is what this double-rounding regression exercises.  Skip if the
	 * environment is not (or cannot report) round-to-nearest -- e.g. arm
	 * under qemu, where fegetround() returns -1 because the FPSCR rounding
	 * field is not read back correctly.
	 */
	if (fegetround() != FE_TONEAREST)
		return REGRESS_SKIP;

	/* The exact reproducer from issue #160. */
	CHECK(fma3(0.9474001f, 4.639901e-7f, -0.24325085f) == -0x1.f22d46p-3f);

	/* A few more triples cross-checked against a hardware FMA. */
	CHECK(fma3(0x1.000002p+0f, 0x1.000002p+0f, -0x1p+0f) == 0x1p-22f);
	CHECK(fma3(0x1.fffffep+0f, 0x1.fffffep+0f, -0x1.fffffcp+1f) == 0x1p-46f);
	CHECK(fma3(0x1.cp+1f, 0x1p+1f, 0x1.8p+0f) == 0x1.1p+3f);
	CHECK(fma3(-0x1.921fb6p+1f, 0x1.5bf0a8p-1f, 0x1p+0f) == -0x1.228b02p+0f);
	CHECK(fma3(0x1.6a09e6p-1f, 0x1.6a09e6p-1f, -0x1p-1f) == -0x1.26055cp-26f);

	return 0;
}
