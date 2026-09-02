/*
 * Regression test for issue #252.  Hard-float RISC-V openlibm used to
 * interpose its FreeBSD-ABI fenv functions on callers using the incompatible
 * Linux ABI.  Use the system <fenv.h> and pull in libopenlibm through
 * isopenlibm().
 */
#include <fenv.h>

#include "regress-util.h"

int isopenlibm(void);

#if !defined(__riscv) || defined(__riscv_float_abi_soft) || \
    !defined(FE_TONEAREST) || !defined(FE_TOWARDZERO) || !defined(FE_UPWARD)

int
main(void)
{
	return REGRESS_SKIP;
}

#else

/*
 * Keep the operands volatile and the division out of line: without
 * -frounding-math the compiler otherwise folds or moves the divisions across
 * the fesetround() calls.
 */
static volatile double one = 1.0, three = 3.0;

static double __attribute__((noinline))
divide(void)
{
	volatile double r = one / three;
	return r;
}

int
main(void)
{
	double nearest, down, up;

	CHECK(isopenlibm() == 1);

	nearest = divide();
	CHECK(fesetround(FE_TOWARDZERO) == 0);
	CHECK(fegetround() == FE_TOWARDZERO);
	down = divide();
	CHECK(fesetround(FE_UPWARD) == 0);
	CHECK(fegetround() == FE_UPWARD);
	up = divide();
	CHECK(fesetround(FE_TONEAREST) == 0);
	CHECK(fegetround() == FE_TONEAREST);

	/* 1/3 is inexact, so the directed modes must bracket the nearest result. */
	CHECK(down < up);
	CHECK(down <= nearest && nearest <= up);
	return 0;
}

#endif
