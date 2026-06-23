/*
 * Regression test for issue #222: powl() is not thread-safe.
 * https://github.com/JuliaMath/openlibm/issues/222
 *
 * The ld80 implementation (ld80/e_powl.c) kept its working scratch values
 * (z, w, W, Wa, Wb, ya, yb, u) in file-scope `static` variables.  When two
 * threads ran powl() at the same time they clobbered each other's scratch,
 * so a call could return arbitrary garbage -- e.g. the reproducer below,
 * powl(0x8.779021e7c2f81b2p+14982L, -0x8.0021b03e1f821c9p-15097L), which is
 * exactly 1, would instead come back as huge values or inf.
 *
 * The fix moves that scratch to local automatic variables.  This test drives
 * several threads, each repeatedly computing a *different* power whose result
 * is known exactly, and fails if any call ever deviates.  With the bug the
 * cross-thread contention is detected essentially every run; after the fix it
 * is rock solid.
 *
 * Skips (exit 77) where it cannot apply: when long double is not the 80-bit
 * extended type this file targets, or when threads cannot be started.
 */
#include <float.h>
#include <math.h>

#include "regress-util.h"

/*
 * ld80-specific: the race was in ld80/e_powl.c, which openlibm builds only on
 * x86.  Restrict to the 80-bit format so every other arch compiles to a bare
 * skip with no reference to powl (it may be absent, breaking the link) and the
 * pthread code is only built where the test actually runs (natively, never
 * under cross-qemu).
 */
#if LDBL_MANT_DIG != 64

int
main(void)
{
	return REGRESS_SKIP;	/* not the 80-bit long double powl */
}

#else

#include <pthread.h>

#define NTHREADS	8
#define NITERS		20000

/*
 * Each job is a (base, exponent, expected) triple chosen so the result is
 * representable exactly in 80-bit long double and goes through the full
 * mantissa/log path of powl() (rather than an early special-case return),
 * which is where the shared scratch lived.  Mixing distinct values across
 * threads is what makes the old race observable.
 */
struct job {
	long double x;
	long double y;
	long double expect;
};

static const struct job jobs[NTHREADS] = {
	/* the exact reproducer from the issue: result is 1 */
	{ 0x8.779021e7c2f81b2p+14982L, -0x8.0021b03e1f821c9p-15097L, 1.0L },
	{ 2.0L,   10.0L, 1024.0L },
	{ 3.0L,    7.0L, 2187.0L },
	{ 1.5L,    4.0L, 5.0625L },
	{ 10.0L,   3.0L, 1000.0L },
	{ 7.0L,    2.0L, 49.0L },
	{ 5.0L,    3.0L, 125.0L },
	{ 2.0L,   -3.0L, 0.125L },
};

static volatile int failed;

static void *
worker(void *arg)
{
	const struct job *j = (const struct job *)arg;
	long i;

	for (i = 0; i < NITERS; i++) {
		long double r = powl(j->x, j->y);
		if (r != j->expect) {
			fprintf(stderr,
			    "  powl(%.20Lg, %.20Lg) = %.20Lg, expected %.20Lg\n",
			    j->x, j->y, r, j->expect);
			failed = 1;
			return NULL;
		}
	}
	return NULL;
}

int
main(void)
{
	pthread_t t[NTHREADS];
	int i, started;

	/* Sanity: every job must be correct single-threaded first. */
	for (i = 0; i < NTHREADS; i++)
		CHECK(powl(jobs[i].x, jobs[i].y) == jobs[i].expect);

	started = 0;
	for (i = 0; i < NTHREADS; i++) {
		if (pthread_create(&t[i], NULL, worker,
		    (void *)&jobs[i]) != 0)
			break;
		started++;
	}

	if (started < 2) {
		/* Could not get real concurrency; nothing meaningful to test. */
		for (i = 0; i < started; i++)
			pthread_join(t[i], NULL);
		return REGRESS_SKIP;
	}

	for (i = 0; i < started; i++)
		pthread_join(t[i], NULL);

	CHECK(failed == 0);
	return 0;
}

#endif /* LDBL_MANT_DIG */
