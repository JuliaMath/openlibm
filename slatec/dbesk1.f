*DECK DBESK1
      DOUBLE PRECISION FUNCTION DBESK1 (X)
C***BEGIN PROLOGUE  DBESK1
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      DOUBLE PRECISION (BESK1-S, DBESK1-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DBESK1(X) calculates the double precision modified (hyperbolic)
C Bessel function of the third kind of order one for double precision
C argument X.  The argument must be large enough that the result does
C not overflow and small enough that the result does not underflow.
C
C Series for BK1        on the interval  0.          to  4.00000E+00
C                                        with weighted error   9.16E-32
C                                         log weighted error  31.04
C                               significant figures required  30.61
C                                    decimal places required  31.64
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DBESI1, DBSK1E, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DBESK1
      DOUBLE PRECISION X, BK1CS(16), XMAX, XMAXT, XMIN, XSML, Y,
     1  D1MACH, DCSEVL, DBESI1, DBSK1E
      LOGICAL FIRST
      SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
      DATA BK1CS(  1) / +.2530022733 8947770532 5311208685 33 D-1     /
      DATA BK1CS(  2) / -.3531559607 7654487566 7238316918 01 D+0     /
      DATA BK1CS(  3) / -.1226111808 2265714823 4790679300 42 D+0     /
      DATA BK1CS(  4) / -.6975723859 6398643501 8129202960 83 D-2     /
      DATA BK1CS(  5) / -.1730288957 5130520630 1765073689 79 D-3     /
      DATA BK1CS(  6) / -.2433406141 5659682349 6007350301 64 D-5     /
      DATA BK1CS(  7) / -.2213387630 7347258558 3152525451 26 D-7     /
      DATA BK1CS(  8) / -.1411488392 6335277610 9583302126 08 D-9     /
      DATA BK1CS(  9) / -.6666901694 1993290060 8537512643 73 D-12    /
      DATA BK1CS( 10) / -.2427449850 5193659339 2631968648 53 D-14    /
      DATA BK1CS( 11) / -.7023863479 3862875971 7837971200 00 D-17    /
      DATA BK1CS( 12) / -.1654327515 5100994675 4910293333 33 D-19    /
      DATA BK1CS( 13) / -.3233834745 9944491991 8933333333 33 D-22    /
      DATA BK1CS( 14) / -.5331275052 9265274999 4666666666 66 D-25    /
      DATA BK1CS( 15) / -.7513040716 2157226666 6666666666 66 D-28    /
      DATA BK1CS( 16) / -.9155085717 6541866666 6666666666 66 D-31    /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DBESK1
      IF (FIRST) THEN
         NTK1 = INITDS (BK1CS, 16, 0.1*REAL(D1MACH(3)))
         XMIN = EXP(MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
         XSML = SQRT(4.0D0*D1MACH(3))
         XMAXT = -LOG(D1MACH(1))
         XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0D0) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X SO SMALL K1 OVERFLOWS', 3, 2)
      Y = 0.D0
      IF (X.GT.XSML) Y = X*X
      DBESK1 = LOG(0.5D0*X)*DBESI1(X) + (0.75D0 + DCSEVL (.5D0*Y-1.D0,
     1  BK1CS, NTK1))/X
      RETURN
C
 20   DBESK1 = 0.D0
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DBESK1',
     +   'X SO BIG K1 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      DBESK1 = EXP(-X) * DBSK1E(X)
C
      RETURN
      END
