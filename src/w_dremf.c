/*
 * dremf() wrapper for remainderf().
 *
 * Written by J.T. Conklin, <jtc@wimsey.com>
 * Placed into the Public Domain, 1994.
 */
/* $FreeBSD: src/lib/msun/src/w_dremf.c,v 1.3 2004/07/28 05:53:18 kan Exp $ */

#include <openlibm_math.h>

#include "math_private.h"

DLLEXPORT float
dremf(float x, float y)
{
	return remainderf(x, y);
}
