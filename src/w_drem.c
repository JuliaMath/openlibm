/*
 * drem() wrapper for remainder().
 *
 * Written by J.T. Conklin, <jtc@wimsey.com>
 * Placed into the Public Domain, 1994.
 */

#include <openlibm_math.h>
#include "math_private.h"

DLLEXPORT double
drem(x, y)
	double x, y;
{
	return remainder(x, y);
}
