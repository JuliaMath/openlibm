#include "cdefs-compat.h"

#include <openlibm.h>

#include "math_private.h"

DLLEXPORT long double
lgammal(long double x)
{
#ifdef OPENLIBM_ONLY_THREAD_SAFE
	int signgam;
#endif

	return (lgammal_r(x, &signgam));
}
