#ifdef __arm__
#include "../arm/fenv.h"
#else
#ifdef __LP64
#include "../amd64/fenv.h"
#else
#include "../i387/fenv.h"
#endif
#endif
