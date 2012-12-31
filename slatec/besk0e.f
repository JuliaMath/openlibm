*DECK BESK0E
      FUNCTION BESK0E (X)
C***BEGIN PROLOGUE  BESK0E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the third kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESK0E-S, DBSK0E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESK0E(X) computes the exponentially scaled modified (hyperbolic)
C Bessel function of third kind of order zero for real argument
C X .GT. 0.0, i.e., EXP(X)*K0(X).
C
C Series for BK0        on the interval  0.          to  4.00000D+00
C                                        with weighted error   3.57E-19
C                                         log weighted error  18.45
C                               significant figures required  17.99
C                                    decimal places required  18.97
C
C Series for AK0        on the interval  1.25000D-01 to  5.00000D-01
C                                        with weighted error   5.34E-17
C                                         log weighted error  16.27
C                               significant figures required  14.92
C                                    decimal places required  16.89
C
C Series for AK02       on the interval  0.          to  1.25000D-01
C                                        with weighted error   2.34E-17
C                                         log weighted error  16.63
C                               significant figures required  14.67
C                                    decimal places required  17.20
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI0, CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESK0E
      DIMENSION BK0CS(11), AK0CS(17), AK02CS(14)
      LOGICAL FIRST
      SAVE BK0CS, AK0CS, AK02CS, NTK0, NTAK0, NTAK02, XSML, FIRST
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
      DATA AK0CS( 1) /   -.0764394790 3327941E0 /
      DATA AK0CS( 2) /   -.0223565260 5699819E0 /
      DATA AK0CS( 3) /    .0007734181 1546938E0 /
      DATA AK0CS( 4) /   -.0000428100 6688886E0 /
      DATA AK0CS( 5) /    .0000030817 0017386E0 /
      DATA AK0CS( 6) /   -.0000002639 3672220E0 /
      DATA AK0CS( 7) /    .0000000256 3713036E0 /
      DATA AK0CS( 8) /   -.0000000027 4270554E0 /
      DATA AK0CS( 9) /    .0000000003 1694296E0 /
      DATA AK0CS(10) /   -.0000000000 3902353E0 /
      DATA AK0CS(11) /    .0000000000 0506804E0 /
      DATA AK0CS(12) /   -.0000000000 0068895E0 /
      DATA AK0CS(13) /    .0000000000 0009744E0 /
      DATA AK0CS(14) /   -.0000000000 0001427E0 /
      DATA AK0CS(15) /    .0000000000 0000215E0 /
      DATA AK0CS(16) /   -.0000000000 0000033E0 /
      DATA AK0CS(17) /    .0000000000 0000005E0 /
      DATA AK02CS( 1) /   -.0120186982 6307592E0 /
      DATA AK02CS( 2) /   -.0091748526 9102569E0 /
      DATA AK02CS( 3) /    .0001444550 9317750E0 /
      DATA AK02CS( 4) /   -.0000040136 1417543E0 /
      DATA AK02CS( 5) /    .0000001567 8318108E0 /
      DATA AK02CS( 6) /   -.0000000077 7011043E0 /
      DATA AK02CS( 7) /    .0000000004 6111825E0 /
      DATA AK02CS( 8) /   -.0000000000 3158592E0 /
      DATA AK02CS( 9) /    .0000000000 0243501E0 /
      DATA AK02CS(10) /   -.0000000000 0020743E0 /
      DATA AK02CS(11) /    .0000000000 0001925E0 /
      DATA AK02CS(12) /   -.0000000000 0000192E0 /
      DATA AK02CS(13) /    .0000000000 0000020E0 /
      DATA AK02CS(14) /   -.0000000000 0000002E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESK0E
      IF (FIRST) THEN
         NTK0 = INITS (BK0CS, 11, 0.1*R1MACH(3))
         NTAK0 = INITS (AK0CS, 17, 0.1*R1MACH(3))
         NTAK02 = INITS (AK02CS, 14, 0.1*R1MACH(3))
         XSML = SQRT (4.0*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.) CALL XERMSG ('SLATEC', 'BESK0E',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.) GO TO 20
C
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESK0E = EXP(X) * (-LOG(0.5*X)*BESI0(X)
     1  - .25 + CSEVL (.5*Y-1., BK0CS, NTK0) )
      RETURN
C
 20   IF (X.LE.8.) BESK0E = (1.25 + CSEVL ((16./X-5.)/3., AK0CS, NTAK0))
     1  / SQRT(X)
      IF (X.GT.8.) BESK0E = (1.25 + CSEVL (16./X-1., AK02CS, NTAK02))
     1  / SQRT(X)
C
      RETURN
      END
