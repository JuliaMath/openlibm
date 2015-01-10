#ifdef OPENLIBM_USE_HOST_FENV_H
#include <fenv.h>
#else /* !OPENLIBM_USE_HOST_FENV_H */

#if defined(__arm__)
#include "../arm/fenv.h"
#elif defined(__x86_64__)
#include "../amd64/fenv.h"
#elif defined(__i386__)
#include "../i387/fenv.h"
#else
#error "Unsupported platform"
#endif

#endif /* OPENLIBM_USE_HOST_FENV_H */
