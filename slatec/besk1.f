*DECK BESK1
      FUNCTION BESK1 (X)
C***BEGIN PROLOGUE  BESK1
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESK1-S, DBESK1-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESK1(X) computes the modified (hyperbolic) Bessel function of third
C kind of order one for real argument X, where X .GT. 0.
C
C Series for BK1        on the interval  0.          to  4.00000D+00
C                                        with weighted error   7.02E-18
C                                         log weighted error  17.15
C                               significant figures required  16.73
C                                    decimal places required  17.67
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI1, BESK1E, CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESK1
      DIMENSION BK1CS(11)
      LOGICAL FIRST
      SAVE BK1CS, NTK1, XMIN, XSML, XMAX, FIRST
      DATA BK1CS( 1) /    .0253002273 389477705E0 /
      DATA BK1CS( 2) /   -.3531559607 76544876E0 /
      DATA BK1CS( 3) /   -.1226111808 22657148E0 /
      DATA BK1CS( 4) /   -.0069757238 596398643E0 /
      DATA BK1CS( 5) /   -.0001730288 957513052E0 /
      DATA BK1CS( 6) /   -.0000024334 061415659E0 /
      DATA BK1CS( 7) /   -.0000000221 338763073E0 /
      DATA BK1CS( 8) /   -.0000000001 411488392E0 /
      DATA BK1CS( 9) /   -.0000000000 006666901E0 /
      DATA BK1CS(10) /   -.0000000000 000024274E0 /
      DATA BK1CS(11) /   -.0000000000 000000070E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESK1
      IF (FIRST) THEN
         NTK1 = INITS (BK1CS, 11, 0.1*R1MACH(3))
         XMIN = EXP (MAX(LOG(R1MACH(1)), -LOG(R1MACH(2))) + .01)
         XSML = SQRT (4.0*R1MACH(3))
         XMAXT = -LOG(R1MACH(1))
         XMAX = XMAXT - 0.5*XMAXT*LOG(XMAXT)/(XMAXT+0.5)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.) CALL XERMSG ('SLATEC', 'BESK1',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'BESK1',
     +   'X SO SMALL K1 OVERFLOWS', 3, 2)
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESK1 = LOG(0.5*X)*BESI1(X) +
     1  (0.75 + CSEVL (.5*Y-1., BK1CS, NTK1))/X
      RETURN
C
 20   BESK1 = 0.
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'BESK1',
     +   'X SO BIG K1 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      BESK1 = EXP(-X) * BESK1E(X)
C
      RETURN
      END
