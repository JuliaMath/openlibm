*DECK FUNDOC
      SUBROUTINE FUNDOC
C***BEGIN PROLOGUE  FUNDOC
C***PURPOSE  Documentation for FNLIB, a collection of routines for
C            evaluating elementary and special functions.
C***LIBRARY   SLATEC
C***CATEGORY  C, Z
C***TYPE      ALL (FUNDOC-A)
C***KEYWORDS  DOCUMENTATION, ELEMENTARY FUNCTIONS, SPECIAL FUNCTIONS
C***AUTHOR  Kahaner, D. K., (NBS)
C***DESCRIPTION
C
C The SLATEC Library --  Elementary And Special Functions
C
C This describes the elementary and special function routines available
C in the SLATEC library.  Most of the these routines were written by
C Wayne Fullerton while at LANL.  Some were written by Don Amos of SNLA.
C There are approximately 63 single precision, 63 double precision and
C 25 complex user callable elementary and special function routines.
C
C The table below gives a breakdown of routines according to their
C function.  Unless otherwise indicated all routines are function
C subprograms.
C                                             Sngl.      Dble.
C Description              Notation           Prec.      Prec.   Complex
C
C         ***Intrinsic Functions and Fundamental Functions***
C Unpack floating point              Call R9UPAK(X,Y,N)  D9UPAK    --
C  number
C Pack floating point                        R9PAK(Y,N)  D9PAK     --
C  number
C Initialize orthogonal               INITS(OS,NOS,ETA)  INITDS    --
C  polynomial series
C Evaluate Chebyshev       summation for  CSEVL(X,CS,N)  DCSEVL    --
C series                  i = 1 to n of
C                          cs(i)*(2*x)**(i-1)
C
C                  ***Elementary Functions***
C Argument = theta in      z = \ z \ *          --         --    CARG(Z)
C  radians                 e**(i * theta)
C Cube root                                   CBRT(X)    DCBRT   CCBRT
C Relative error exponen-  ((e**x) -1) / x    EXPREL(X)  DEXPRL  CEXPRL
C  tial from first order
C Common logarithm         log to the base 10   --         --  CLOG10(Z)
C                          of z
C Relative error logarithm ln(1 + x)          ALNREL(X)  DLNREL  CLNREL
C Relative error logarithm (ln(1 + x) - x     R9LN2R(X)  D9LN2R  C9LN2R
C from second order        + x**2/2) / x**3
C               ***Trigonometric and Hyperbolic Functions***
C Tangent                  tan z                --         --    CTAN(Z)
C Cotangent                cot x              COT(X)     DCOT    CCOT
C Sine x in degrees        sin((2*pi*x)/360)  SINDG(X)   DSINDG    --
C Cosine x in degrees      cos((2*pi*x)/360)  COSDG(X)   DCOSDG    --
C Arc sine                 arcsin (z)           --         --   CASIN(Z)
C Arc cosine               arccos (z)           --         --   CACOS(Z)
C Arc tangent              arctan (z)           --         --   CATAN(Z)
C Quadrant correct         arctan (z1/z2)       --         -- CATAN2(Z1,
C  arc tangent                                                       Z2)
C Hyperbolic sine          sinh z               --         --   CSINH(Z)
C Hyperbolic cosine        cosh z               --         --   CCOSH(Z)
C Hyperbolic tangent       tanh z               --         --   CTANH(Z)
C Arc hyperbolic sine      arcsinh (x)        ASINH(X)   DASINH  CASINH
C Arc hyperbolic cosine    arccosh (x)        ACOSH(X)   DACOSH  CACOSH
C Arc hyperbolic tangent   arctanh (x)        ATANH(X)   DATANH  CATANH
C Relative error arc       (arctan (x) - x)   R9ATN1(X)  D9ATN1    --
C  tangent from first order   / x**3
C              ***Exponential Integrals and Related Functions***
C Exponential integral     Ei(x) = (minus)    EI(X)      DEI       --
C                          the integral from
C                          -x to infinity of
C                            (e**-t / t)dt
C Exponential integral     E sub 1 (x) =      E1(X)      DE1       --
C                          the integral from x
C                            to infinity of
C                          (e**-t / t) dt
C Logarithmic integral     li(x) = the        ALI(X)     DLI       --
C                          integral from 0 to
C                          x of (1 / ln t) dt
C   Sequences of exponential integrals.
C   M values are computed where
C   k=0,1,...M-1 and n>=1
C Exponential integral     E sub n+k (x) Call EXINT(X,   DEXINT    --
C                        =the integral from   N,KODE,M,TOL,
C                         1 to infinity of    EN,IERR)
C                       (e**(-x*t)/t**(n+k))dt
C                 ***Gamma Functions and Related Functions***
C Factorial                n!                 FAC(N)     DFAC      --
C Binomial                 n!/(m!*(n-m)!)     BINOM(N,M) DBINOM    --
C Gamma                    gamma(x)           GAMMA(X)   DGAMMA  CGAMMA
C Gamma(x) under and                     Call GAMLIM(    DGAMLM    --
C  overflow limits                           XMIN,XMAX)
C Reciprocal gamma         1 / gamma(x)       GAMR(X)    DGAMR   CGAMR
C Log abs gamma            ln \gamma(x)\      ALNGAM(X)  DLNGAM    --
C Log gamma                ln gamma(z)          --         --    CLNGAM
C Log abs gamma       g = ln \gamma(x)\  Call ALGAMS(X,  DLGAMS    --
C with sign           s = sign gamma(x)      G,S)
C Incomplete gamma         gamma(a,x) =       GAMI(A,X)  DGAMI     --
C                          the integral from
C                          0 to x of
C                         (t**(a-1) * e**-t)dt
C Complementary            gamma(a,x) =       GAMIC(A,X) DGAMIC    --
C  incomplete gamma        the integral from
C                          x to infinity of
C                         (t**(a-1) * e**-t)dt
C Tricomi's             gamma super star(a,x) GAMIT(A,X) DGAMIT    --
C  incomplete gamma        = x**-a *
C                         incomplete gamma(a,x)
C                          / gamma(a)
C Psi (Digamma)            psi(x) = gamma'(x) PSI(X)     DPSI    CPSI
C                          / gamma(x)
C Pochhammer's         (a) sub x = gamma(a+x) POCH(A,X)  DPOCH     --
C  generalized symbol      / gamma(a)
C Pochhammer's symbol    ((a) sub x -1) / x   POCH1(A,X) DPOCH1    --
C  from first order
C Beta                     b(a,b) = (gamma(a) BETA(A,B)  DBETA   CBETA
C                          * gamma(b))
C                          / gamma(a+b)
C                           = the integral
C                           from 0 to 1 of
C                           (t**(a-1) *
C                           (1-t)**(b-1))dt
C Log beta                 ln b(a,b)         ALBETA(A,B) DLBETA  CLBETA
C Incomplete beta          i sub x (a,b) =  BETAI(X,A,B) DBETAI    __
C                          b sub x (a,b) / b(a,b)
C                           = 1 / b(a,b) *
C                          the integral
C                          from 0 to x of
C                          (t**(a-1) *
C                          (1-t)**(b-1))dt
C Log gamma correction     ln gamma(x) -      R9LGMC(X)  D9LGMC  C9LGMC
C  term when Stirling's    (ln(2 * pi))/2 -
C  approximation is valid  (x - 1/2) * ln(x) + x
C                ***Error Functions and Fresnel Integrals***
C Error function           erf x = (2 /       ERF(X)     DERF      --
C                          square root of pi) *
C                          the integral from
C                          0 to x of
C                          e**(-t**2)dt
C Complementary            erfc x = (2 /      ERFC(X)    DERFC     --
C  error function          square root of pi) *
C                          the integral from
C                          x to infinity of
C                          e**(-t**2)dt
C Dawson's function        F(x) = e**(-x**2)  DAWS(X)    DDAWS     --
C                          * the integral from
C                          from 0 to x of
C                          e**(t**2)dt
C                         ***Bessel Functions***
C   Bessel functions of special integer order
C First kind, order zero   J sub 0 (x)        BESJ0(X)   DBESJ0    --
C First kind, order one    J sub 1 (x)        BESJ1(X)   DBESJ1    --
C Second kind, order zero  Y sub 0 (x)        BESY0(X)   DBESY0    --
C Second kind, order one   Y sub 1 (x)        BESY1(X)   DBESY1    --
C   Modified (hyperbolic) Bessel functions of special integer order
C First kind, order zero   I sub 0 (x)        BESI0(X)   DBESI0    --
C First kind, order one    I sub 1 (x)        BESI1(X)   DBESI1    --
C Third kind, order zero   K sub 0 (x)        BESK0(X)   DBESK0    --
C Third kind, order one    K sub 1 (x)        BESK1(X)   DBESK1    --
C   Modified (hyperbolic) Bessel functions of special integer order
C   scaled by an exponential
C First kind, order zero   e**-\x\ * I sub 0(x) BESI0E(X) DBSI0E   --
C First kind, order one    e**-\x\ * I sub 1(x) BESI1E(X) DBSI1E   --
C Third kind, order zero   e**x * K sub 0 (x)   BESK0E(X) DBSK0E   --
C Third kind, order one    e**x * K sub 1 (x)   BESK1E(X) DBSK1E   --
C   Sequences of Bessel functions of general order.
C   N values are computed where  k = 1,2,...N and v .ge. 0.
C Modified first kind      I sub v+k-1 (x) Call BESI(X,   DBESI    --
C                          optional scaling  ALPHA,KODE,N,
C                          by e**(-x)        Y,NZ)
C First kind               J sub v+k-1 (x) Call BESJ(X,   DBESJ    --
C                                            ALPHA,N,Y,NZ)
C Second kind              Y sub v+k-1 (x) Call BESY(X,   DBESY    --
C                                            FNU,N,Y)
C Modified third kind      K sub v+k-1 (x) Call BESK(X,   DBESK    --
C                          optional scaling  FNU,KODE,N,Y,
C                          by e**(x)         NZ)
C   Sequences of Bessel functions.  \N\ values are computed where
C   I = 0, 1, 2, ..., N-1  for N > 0  or I = 0, -1, -2, ..., N+1
C   for N < 0.
C Modified third kind      K sub v+i (x)   Call BESKS(    DBESKS   --
C                                           XNU,X,N,BK)
C   Sequences of Bessel functions scaled by an exponential.
C   \N\ values are computed where  I = 0, 1, 2, ..., N-1
C   for N > 0  or  I = 0, -1, -2, ..., N+1  for N < 0.
C Modified third kind      e**x *         Call BESKES(    DBSKES   --
C                          K sub v+i (x)     XNU,X,N,BK)
C                ***Bessel Functions of Fractional Order***
C   Airy functions
C Airy                     Ai(x)              AI(X)      DAI       --
C Bairy                    Bi(x)              BI(X)      DBI       --
C   Exponentially scaled Airy functions
C Airy                     Ai(x), x <= 0      AIE(X)     DAIE      --
C                          exp(2/3 * x**(3/2))
C                          * Ai(x), x >= 0
C Bairy                    Bi(x), x <= 0      BIE(X)     DBIE      --
C                          exp(-2/3 * x**(3/2))
C                          * Bi(x), x >= 0
C                 ***Confluent Hypergeometric Functions***
C Confluent                U(a,b,x)           CHU(A,B,X) DCHU      --
C  hypergeometric
C                     ***Miscellaneous Functions***
C Spence                   s(x) = - the       SPENC(X)   DSPENC    --
C  dilogarithm             integral from
C                          0 to x of
C                          ((ln \1-y\) / y)dy
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801015  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Routine name changed from FNLIBD to FUNDOC.  (WRB)
C   900723  PURPOSE section revised.  (WRB)
C***END PROLOGUE  FUNDOC
C***FIRST EXECUTABLE STATEMENT  FUNDOC
      RETURN
      END
