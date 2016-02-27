#ifndef _OPENLIBM_FENV_H

#define _OPENLIBM_FENV_H

#ifdef OPENLIBM_USE_HOST_FENV_H
#include <fenv.h>
#else /* !OPENLIBM_USE_HOST_FENV_H */

#if defined(__arm__)
#include <openlibm_fenv_arm.h>
#elif defined(__aarch64__)
#include <openlibm_fenv_aarch64.h>
#elif defined(__i386__) || defined(__x86_64__)
#include <openlibm_fenv_x86.h>
#elif defined(__powerpc__)
#include <openlibm_fenv_powerpc.h>
#else
#error "Unsupported platform: Try -DOPENLIBM_USE_HOST_FENV_H"
#endif

#endif /* OPENLIBM_USE_HOST_FENV_H */

#endif /* _OPENLIBM_FENV_H */
