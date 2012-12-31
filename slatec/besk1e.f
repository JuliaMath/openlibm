*DECK BESK1E
      FUNCTION BESK1E (X)
C***BEGIN PROLOGUE  BESK1E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the third kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESK1E-S, DBSK1E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
C             THIRD KIND
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESK1E(X) computes the exponentially scaled modified (hyperbolic)
C Bessel function of third kind of order one for real argument
C X .GT. 0.0, i.e., EXP(X)*K1(X).
C
C Series for BK1        on the interval  0.          to  4.00000D+00
C                                        with weighted error   7.02E-18
C                                         log weighted error  17.15
C                               significant figures required  16.73
C                                    decimal places required  17.67
C
C Series for AK1        on the interval  1.25000D-01 to  5.00000D-01
C                                        with weighted error   6.06E-17
C                                         log weighted error  16.22
C                               significant figures required  15.41
C                                    decimal places required  16.83
C
C Series for AK12       on the interval  0.          to  1.25000D-01
C                                        with weighted error   2.58E-17
C                                         log weighted error  16.59
C                               significant figures required  15.22
C                                    decimal places required  17.16
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI1, CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESK1E
      DIMENSION BK1CS(11), AK1CS(17), AK12CS(14)
      LOGICAL FIRST
      SAVE BK1CS, AK1CS, AK12CS, NTK1, NTAK1, NTAK12, XMIN, XSML,
     1 FIRST
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
      DATA AK1CS( 1) /    .2744313406 973883E0 /
      DATA AK1CS( 2) /    .0757198995 3199368E0 /
      DATA AK1CS( 3) /   -.0014410515 5647540E0 /
      DATA AK1CS( 4) /    .0000665011 6955125E0 /
      DATA AK1CS( 5) /   -.0000043699 8470952E0 /
      DATA AK1CS( 6) /    .0000003540 2774997E0 /
      DATA AK1CS( 7) /   -.0000000331 1163779E0 /
      DATA AK1CS( 8) /    .0000000034 4597758E0 /
      DATA AK1CS( 9) /   -.0000000003 8989323E0 /
      DATA AK1CS(10) /    .0000000000 4720819E0 /
      DATA AK1CS(11) /   -.0000000000 0604783E0 /
      DATA AK1CS(12) /    .0000000000 0081284E0 /
      DATA AK1CS(13) /   -.0000000000 0011386E0 /
      DATA AK1CS(14) /    .0000000000 0001654E0 /
      DATA AK1CS(15) /   -.0000000000 0000248E0 /
      DATA AK1CS(16) /    .0000000000 0000038E0 /
      DATA AK1CS(17) /   -.0000000000 0000006E0 /
      DATA AK12CS( 1) /    .0637930834 3739001E0 /
      DATA AK12CS( 2) /    .0283288781 3049721E0 /
      DATA AK12CS( 3) /   -.0002475370 6739052E0 /
      DATA AK12CS( 4) /    .0000057719 7245160E0 /
      DATA AK12CS( 5) /   -.0000002068 9392195E0 /
      DATA AK12CS( 6) /    .0000000097 3998344E0 /
      DATA AK12CS( 7) /   -.0000000005 5853361E0 /
      DATA AK12CS( 8) /    .0000000000 3732996E0 /
      DATA AK12CS( 9) /   -.0000000000 0282505E0 /
      DATA AK12CS(10) /    .0000000000 0023720E0 /
      DATA AK12CS(11) /   -.0000000000 0002176E0 /
      DATA AK12CS(12) /    .0000000000 0000215E0 /
      DATA AK12CS(13) /   -.0000000000 0000022E0 /
      DATA AK12CS(14) /    .0000000000 0000002E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESK1E
      IF (FIRST) THEN
         NTK1 = INITS (BK1CS, 11, 0.1*R1MACH(3))
         NTAK1 = INITS (AK1CS, 17, 0.1*R1MACH(3))
         NTAK12 = INITS (AK12CS, 14, 0.1*R1MACH(3))
C
         XMIN = EXP (MAX(LOG(R1MACH(1)), -LOG(R1MACH(2))) + .01)
         XSML = SQRT (4.0*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.) CALL XERMSG ('SLATEC', 'BESK1E',
     +   'X IS ZERO OR NEGATIVE', 2, 2)
      IF (X.GT.2.0) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'BESK1E',
     +   'X SO SMALL K1 OVERFLOWS', 3, 2)
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESK1E = EXP(X) * (LOG(0.5*X)*BESI1(X) +
     1  (0.75 + CSEVL (.5*Y-1., BK1CS, NTK1))/X )
      RETURN
C
 20   IF (X.LE.8.) BESK1E = (1.25 + CSEVL ((16./X-5.)/3., AK1CS, NTAK1))
     1  / SQRT(X)
      IF (X.GT.8.) BESK1E = (1.25 + CSEVL (16./X-1., AK12CS, NTAK12))
     1  / SQRT(X)
C
      RETURN
      END
