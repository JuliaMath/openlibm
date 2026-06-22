/*
 * Regression test for issue #211: loss of precision in powf.
 * https://github.com/JuliaMath/openlibm/issues/211
 *
 * Fixed in 98f8713 by importing the latest msun powf; this pins the exact
 * result so the regression cannot silently return.
 */
#include <math.h>

#include "regress-util.h"

int
main(void)
{
	float x = 0xd.65874p-4f;
	float y = 4.0f;
	float z = powf(x, y);

	CHECK(z == 0x1.f74424p-2f);
	return 0;
}
