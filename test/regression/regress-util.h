/*
 * Minimal helpers for openlibm standalone regression tests.
 *
 * Each regression test is a self-contained program named after the issue it
 * covers (test/regression/test-<issue>.c).  The harness in test/Makefile builds
 * and runs every file in this directory and interprets the exit status as:
 *
 *     0   the test passed
 *     77  the test does not apply on this platform and was skipped
 *     *   the test failed
 *
 * Use REGRESS_SKIP for the skip status and CHECK() for assertions so a failing
 * test prints which condition broke and where.
 */
#ifndef OPENLIBM_REGRESS_UTIL_H
#define OPENLIBM_REGRESS_UTIL_H

#include <stdio.h>

#define REGRESS_SKIP 77

#define CHECK(cond)							\
	do {								\
		if (!(cond)) {						\
			fprintf(stderr, "  check failed: %s (%s:%d)\n",	\
			    #cond, __FILE__, __LINE__);			\
			return 1;					\
		}							\
	} while (0)

#endif /* OPENLIBM_REGRESS_UTIL_H */
