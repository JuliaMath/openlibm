include ../Make.inc

SRCS= \
	e_acos.o e_acosf.o e_acosh.o e_acoshf.o e_asin.o e_asinf.o \
	e_atan2.o e_atan2f.o e_atanh.o e_atanhf.o e_cosh.o e_coshf.o e_exp.o \
	e_expf.o e_fmod.o e_fmodf.o e_gamma.o e_gamma_r.o e_gammaf.o \
	e_gammaf_r.o e_hypot.o e_hypotf.o e_j0.o e_j0f.o e_j1.o e_j1f.o \
	e_jn.o e_jnf.o e_lgamma.o e_lgamma_r.o e_lgammaf.o e_lgammaf_r.o \
	e_log.o e_log10.o e_log10f.o e_log2.o e_log2f.o e_logf.o \
	e_pow.o e_powf.o e_remainder.o e_remainderf.o e_scalb.o e_scalbf.o \
	e_rem_pio2.o e_rem_pio2f.o \
	e_sinh.o e_sinhf.o e_sqrt.o e_sqrtf.o \
	k_cos.o k_exp.o k_expf.o k_rem_pio2.o k_sin.o k_tan.o \
	k_cosf.o k_sinf.o k_tanf.o \
	s_asinh.o s_asinhf.o s_atan.o s_atanf.o s_carg.o s_cargf.o s_cargl.o \
	s_cbrt.o s_cbrtf.o s_ceil.o s_ceilf.o \
	s_copysign.o s_copysignf.o s_cos.o s_cosf.o \
	s_csqrt.o s_csqrtf.o s_erf.o s_erff.o \
	s_exp2.o s_exp2f.o s_expm1.o s_expm1f.o s_fabs.o s_fabsf.o s_fdim.o \
	s_finite.o s_finitef.o \
	s_floor.o s_floorf.o s_fma.o s_fmaf.o \
	s_fmax.o s_fmaxf.o s_fmaxl.o s_fmin.o \
	s_fminf.o s_fminl.o s_fpclassify.o \
	s_frexp.o s_frexpf.o s_ilogb.o s_ilogbf.o \
	s_ilogbl.o s_isinf.o s_isfinite.o s_isnormal.o s_isnan.o \
	s_llrint.o s_llrintf.o s_llround.o s_llroundf.o s_llroundl.o \
	s_log1p.o s_log1pf.o s_logb.o s_logbf.o s_lrint.o s_lrintf.o \
	s_lround.o s_lroundf.o s_lroundl.o s_modf.o s_modff.o \
	s_nan.o s_nearbyint.o s_nextafter.o s_nextafterf.o \
	s_nexttowardf.o s_remquo.o s_remquof.o \
	s_rint.o s_rintf.o s_round.o s_roundf.o s_roundl.o \
	s_scalbln.o s_scalbn.o s_scalbnf.o s_signbit.o \
	s_signgam.o s_significand.o s_significandf.o s_sin.o s_sinf.o \
	s_tan.o s_tanf.o s_tanh.o s_tanhf.o s_tgammaf.o s_trunc.o s_truncf.o \
	w_cabs.o w_cabsf.o w_drem.o w_dremf.o

# C99 long double functions
SRCS+=	s_copysignl.o s_fabsl.o s_llrintl.o s_lrintl.o s_modfl.o

# If long double != double use these; otherwise, we alias the double versions.
SRCS+=	e_acosl.o e_asinl.o e_atan2l.o e_fmodl.o \
	e_hypotl.o e_remainderl.o e_sqrtl.o \
	s_atanl.o s_ceill.o s_cosl.o s_cprojl.o \
	s_csqrtl.o s_floorl.o s_fmal.o \
	s_frexpl.o s_logbl.o s_nexttoward.o \
	s_remquol.o \
	s_sinl.o s_tanl.o s_truncl.o w_cabsl.o \
	s_nextafterl.o s_rintl.o s_scalbnl.o 
#	s_cbrtl.o

# C99 complex functions
SRCS+=	s_ccosh.o s_ccoshf.o s_cexp.o s_cexpf.o \
	s_cimag.o s_cimagf.o s_cimagl.o \
	s_conj.o s_conjf.o s_conjl.o \
	s_cproj.o s_cprojf.o s_creal.o s_crealf.o s_creall.o \
	s_csinh.o s_csinhf.o s_ctanh.o s_ctanhf.o

all: $(SRCS)
