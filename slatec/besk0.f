*DECK BESK0
      FUNCTION BESK0 (X)
C***BEGIN PROLOGUE  BESK0
C***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
C            third kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESK0-S, DBESK0-D)
C***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESK0(X) calculates the modified (hyperbolic) Bessel function
C of the third kind of order zero for real argument X .GT. 0.0.
C
C Series for BK0        on the interval  0.          to  4.00000D+00
C                                        with weighted error   3.57E-19
C                                         log weighted error  18.45
C                               significant figures required  17.99
C                                    decimal places required  18.97
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI0, BESK0E, CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESK0
      DIMENSION BK0CS(11)
      LOGICAL FIRST
      SAVE BK0CS, NTK0, XSML, XMAX, FIRST
      DATA BK0CS( 1) /   -.0353273932 3390276872E0 /
      DATA BK0CS( 2) /    .3442898999 246284869E0 /
      DATA BK0CS( 3) /    .0359799365 1536150163E0 /
      DATA BK0CS( 4) /    .0012646154 1144692592E0 /
      DATA BK0CS( 5) /    .0000228621 2103119451E0 /
      DATA BK0CS( 6) /    .0000002534 7910790261E0 /
      DATA BK0CS( 7) /    .0000000019 0451637722E0 /
      DATA BK0CS( 8) /    .0000000000 1034969525E0 /
      DATA BK0CS( 9) /    .0000000000 0004259816E0 /
      DATA BK0CS(10) /    .0000000000 0000013744E0 /
      DATA BK0CS(11) /    .0000000000 0000000035E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESK0
      IF (FIRST) THEN
         NTK0 = INITS (BK0CS, 11, 0.1*R1MACH(3))
         XSML = SQRT (4.0*R1MACH(3))
         XMAXT = -LOG(R1MACH(1))
         XMAX = XMAXT - 0.5*XMAXT*LOG(XMAXT)/(XMAXT+0.5) - 0.01
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.) CALL XERMSG ('SLATEC', 'BESK0',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.) GO TO 20
C
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESK0 = -LOG(0.5*X)*BESI0(X) - .25 + CSEVL (.5*Y-1., BK0CS, NTK0)
      RETURN
C
 20   BESK0 = 0.
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'BESK0',
     +   'X SO BIG K0 UNDERFLOWS', 1, 1)
      IF (X.GT.XMAX) RETURN
C
      BESK0 = EXP(-X) * BESK0E(X)
C
      RETURN
      END
