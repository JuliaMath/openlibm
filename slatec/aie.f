*DECK AIE
      FUNCTION AIE (X)
C***BEGIN PROLOGUE  AIE
C***PURPOSE  Calculate the Airy function for a negative argument and an
C            exponentially scaled Airy function for a non-negative
C            argument.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10D
C***TYPE      SINGLE PRECISION (AIE-S, DAIE-D)
C***KEYWORDS  EXPONENTIALLY SCALED AIRY FUNCTION, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C AIE(X) computes the exponentially scaled Airy function for
C non-negative X.  It evaluates AI(X) for X .LE. 0.0 and
C EXP(ZETA)*AI(X) for X .GE. 0.0 where ZETA = (2.0/3.0)*(X**1.5).
C
C Series for AIF        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.09E-19
C                                         log weighted error  18.96
C                               significant figures required  17.76
C                                    decimal places required  19.44
C
C Series for AIG        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.51E-17
C                                         log weighted error  16.82
C                               significant figures required  15.19
C                                    decimal places required  17.27
C
C Series for AIP        on the interval  0.          to  1.00000D+00
C                                        with weighted error   5.10E-17
C                                         log weighted error  16.29
C                               significant figures required  14.41
C                                    decimal places required  17.06
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  AIE
      DIMENSION AIFCS(9), AIGCS(8), AIPCS(34)
      LOGICAL FIRST
      SAVE AIFCS, AIGCS, AIPCS, NAIF, NAIG,
     1 NAIP, X3SML, X32SML, XBIG, FIRST
      DATA AIFCS( 1) /   -.0379713584 9666999750E0 /
      DATA AIFCS( 2) /    .0591918885 3726363857E0 /
      DATA AIFCS( 3) /    .0009862928 0577279975E0 /
      DATA AIFCS( 4) /    .0000068488 4381907656E0 /
      DATA AIFCS( 5) /    .0000000259 4202596219E0 /
      DATA AIFCS( 6) /    .0000000000 6176612774E0 /
      DATA AIFCS( 7) /    .0000000000 0010092454E0 /
      DATA AIFCS( 8) /    .0000000000 0000012014E0 /
      DATA AIFCS( 9) /    .0000000000 0000000010E0 /
      DATA AIGCS( 1) /    .0181523655 8116127E0 /
      DATA AIGCS( 2) /    .0215725631 6601076E0 /
      DATA AIGCS( 3) /    .0002567835 6987483E0 /
      DATA AIGCS( 4) /    .0000014265 2141197E0 /
      DATA AIGCS( 5) /    .0000000045 7211492E0 /
      DATA AIGCS( 6) /    .0000000000 0952517E0 /
      DATA AIGCS( 7) /    .0000000000 0001392E0 /
      DATA AIGCS( 8) /    .0000000000 0000001E0 /
      DATA AIPCS( 1) /   -.0187519297 793868E0 /
      DATA AIPCS( 2) /   -.0091443848 250055E0 /
      DATA AIPCS( 3) /    .0009010457 337825E0 /
      DATA AIPCS( 4) /   -.0001394184 127221E0 /
      DATA AIPCS( 5) /    .0000273815 815785E0 /
      DATA AIPCS( 6) /   -.0000062750 421119E0 /
      DATA AIPCS( 7) /    .0000016064 844184E0 /
      DATA AIPCS( 8) /   -.0000004476 392158E0 /
      DATA AIPCS( 9) /    .0000001334 635874E0 /
      DATA AIPCS(10) /   -.0000000420 735334E0 /
      DATA AIPCS(11) /    .0000000139 021990E0 /
      DATA AIPCS(12) /   -.0000000047 831848E0 /
      DATA AIPCS(13) /    .0000000017 047897E0 /
      DATA AIPCS(14) /   -.0000000006 268389E0 /
      DATA AIPCS(15) /    .0000000002 369824E0 /
      DATA AIPCS(16) /   -.0000000000 918641E0 /
      DATA AIPCS(17) /    .0000000000 364278E0 /
      DATA AIPCS(18) /   -.0000000000 147475E0 /
      DATA AIPCS(19) /    .0000000000 060851E0 /
      DATA AIPCS(20) /   -.0000000000 025552E0 /
      DATA AIPCS(21) /    .0000000000 010906E0 /
      DATA AIPCS(22) /   -.0000000000 004725E0 /
      DATA AIPCS(23) /    .0000000000 002076E0 /
      DATA AIPCS(24) /   -.0000000000 000924E0 /
      DATA AIPCS(25) /    .0000000000 000417E0 /
      DATA AIPCS(26) /   -.0000000000 000190E0 /
      DATA AIPCS(27) /    .0000000000 000087E0 /
      DATA AIPCS(28) /   -.0000000000 000040E0 /
      DATA AIPCS(29) /    .0000000000 000019E0 /
      DATA AIPCS(30) /   -.0000000000 000009E0 /
      DATA AIPCS(31) /    .0000000000 000004E0 /
      DATA AIPCS(32) /   -.0000000000 000002E0 /
      DATA AIPCS(33) /    .0000000000 000001E0 /
      DATA AIPCS(34) /   -.0000000000 000000E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  AIE
      IF (FIRST) THEN
         ETA = 0.1*R1MACH(3)
         NAIF  = INITS (AIFCS, 9, ETA)
         NAIG  = INITS (AIGCS, 8, ETA)
         NAIP  = INITS (AIPCS, 34, ETA)
C
         X3SML = ETA**0.3333
         X32SML = 1.3104*X3SML**2
         XBIG = R1MACH(2)**0.6666
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GE.(-1.0)) GO TO 20
      CALL R9AIMP (X, XM, THETA)
      AIE = XM * COS(THETA)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 30
      Z = 0.0
      IF (ABS(X).GT.X3SML) Z = X**3
      AIE = 0.375 + (CSEVL (Z, AIFCS, NAIF) - X*(0.25 +
     1  CSEVL (Z, AIGCS, NAIG)) )
      IF (X.GT.X32SML) AIE = AIE * EXP(2.0*X*SQRT(X)/3.0)
      RETURN
C
 30   SQRTX = SQRT(X)
      Z = -1.0
      IF (X.LT.XBIG) Z = 2.0/(X*SQRTX) - 1.0
      AIE = (.28125 + CSEVL (Z, AIPCS, NAIP))/SQRT(SQRTX)
      RETURN
C
      END
