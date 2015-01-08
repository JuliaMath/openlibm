#include "cdefs-compat.h"
#include "openlibm.h"
#include "math_private.h"

extern int signgam;

DLLEXPORT long double
lgammal(long double x)
{

	return (lgammal_r(x, &signgam));
}
