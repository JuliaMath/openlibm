
/* @(#)e_lgamma.c 1.3 95/01/18 */
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice 
 * is preserved.
 * ====================================================
 *
 */

#include "cdefs-compat.h"
//__FBSDID("$FreeBSD: src/lib/msun/src/e_lgamma.c,v 1.9 2008/02/22 02:30:35 das Exp $");

/* lgamma(x)
 * Return the logarithm of the Gamma function of x.
 *
 * Method: call lgamma_r
 */

#include <openlibm_math.h>

#include "math_private.h"

OLM_DLLEXPORT double
lgamma(double x)
{
#ifdef OPENLIBM_ONLY_THREAD_SAFE
	int signgam;
#endif

	return lgamma_r(x,&signgam);
}
