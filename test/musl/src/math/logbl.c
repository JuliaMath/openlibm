#include <stdint.h>
#include <stdio.h>
#include "mtest.h"

static struct l_l t[] = {
#if LDBL_MANT_DIG == 53
#include "sanity/logb.h"
#include "special/logb.h"

#elif LDBL_MANT_DIG == 64
#include "sanity/logbl.h"
#include "special/logbl.h"

#endif
};

int main(void)
{
	#pragma STDC FENV_ACCESS ON
	long double y;
	float d;
	int e, i, err = 0;
	struct l_l *p;

	for (i = 0; i < sizeof t/sizeof *t; i++) {
		p = t + i;

		if (p->r < 0)
			continue;
		fesetround(p->r);
		feclearexcept(FE_ALL_EXCEPT);
		y = logbl(p->x);
		e = fetestexcept(INEXACT|INVALID|DIVBYZERO|UNDERFLOW|OVERFLOW);

		if (!checkexceptall(e, p->e, p->r)) {
			printf("%s:%d: bad fp exception: %s logbl(%La)=%La, want %s",
				p->file, p->line, rstr(p->r), p->x, p->y, estr(p->e));
			printf(" got %s\n", estr(e));
			err++;
		}
		d = ulperrl(y, p->y, p->dy);
		if (!checkcr(y, p->y, p->r)) {
			printf("%s:%d: %s logbl(%La) want %La got %La ulperr %.3f = %a + %a\n",
				p->file, p->line, rstr(p->r), p->x, p->y, y, d, d-p->dy, p->dy);
			err++;
		}
	}
	return !!err;
}
