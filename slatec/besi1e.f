*DECK BESI1E
      FUNCTION BESI1E (X)
C***BEGIN PROLOGUE  BESI1E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESI1E-S, DBSI1E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
C             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
C             ORDER ONE, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESI1E(X) calculates the exponentially scaled modified (hyperbolic)
C Bessel function of the first kind of order one for real argument X;
C i.e., EXP(-ABS(X))*I1(X).
C
C Series for BI1        on the interval  0.          to  9.00000D+00
C                                        with weighted error   2.40E-17
C                                         log weighted error  16.62
C                               significant figures required  16.23
C                                    decimal places required  17.14
C
C Series for AI1        on the interval  1.25000D-01 to  3.33333D-01
C                                        with weighted error   6.98E-17
C                                         log weighted error  16.16
C                               significant figures required  14.53
C                                    decimal places required  16.82
C
C Series for AI12       on the interval  0.          to  1.25000D-01
C                                        with weighted error   3.55E-17
C                                         log weighted error  16.45
C                               significant figures required  14.69
C                                    decimal places required  17.12
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890210  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  BESI1E
      DIMENSION BI1CS(11), AI1CS(21), AI12CS(22)
      LOGICAL FIRST
      SAVE BI1CS, AI1CS, AI12CS, NTI1, NTAI1, NTAI12, XMIN, XSML, FIRST
      DATA BI1CS( 1) /   -.0019717132 61099859E0 /
      DATA BI1CS( 2) /    .4073488766 7546481E0 /
      DATA BI1CS( 3) /    .0348389942 99959456E0 /
      DATA BI1CS( 4) /    .0015453945 56300123E0 /
      DATA BI1CS( 5) /    .0000418885 21098377E0 /
      DATA BI1CS( 6) /    .0000007649 02676483E0 /
      DATA BI1CS( 7) /    .0000000100 42493924E0 /
      DATA BI1CS( 8) /    .0000000000 99322077E0 /
      DATA BI1CS( 9) /    .0000000000 00766380E0 /
      DATA BI1CS(10) /    .0000000000 00004741E0 /
      DATA BI1CS(11) /    .0000000000 00000024E0 /
      DATA AI1CS( 1) /   -.0284674418 1881479E0 /
      DATA AI1CS( 2) /   -.0192295323 1443221E0 /
      DATA AI1CS( 3) /   -.0006115185 8579437E0 /
      DATA AI1CS( 4) /   -.0000206997 1253350E0 /
      DATA AI1CS( 5) /    .0000085856 1914581E0 /
      DATA AI1CS( 6) /    .0000010494 9824671E0 /
      DATA AI1CS( 7) /   -.0000002918 3389184E0 /
      DATA AI1CS( 8) /   -.0000000155 9378146E0 /
      DATA AI1CS( 9) /    .0000000131 8012367E0 /
      DATA AI1CS(10) /   -.0000000014 4842341E0 /
      DATA AI1CS(11) /   -.0000000002 9085122E0 /
      DATA AI1CS(12) /    .0000000001 2663889E0 /
      DATA AI1CS(13) /   -.0000000000 1664947E0 /
      DATA AI1CS(14) /   -.0000000000 0166665E0 /
      DATA AI1CS(15) /    .0000000000 0124260E0 /
      DATA AI1CS(16) /   -.0000000000 0027315E0 /
      DATA AI1CS(17) /    .0000000000 0002023E0 /
      DATA AI1CS(18) /    .0000000000 0000730E0 /
      DATA AI1CS(19) /   -.0000000000 0000333E0 /
      DATA AI1CS(20) /    .0000000000 0000071E0 /
      DATA AI1CS(21) /   -.0000000000 0000006E0 /
      DATA AI12CS( 1) /    .0285762350 1828014E0 /
      DATA AI12CS( 2) /   -.0097610974 9136147E0 /
      DATA AI12CS( 3) /   -.0001105889 3876263E0 /
      DATA AI12CS( 4) /   -.0000038825 6480887E0 /
      DATA AI12CS( 5) /   -.0000002512 2362377E0 /
      DATA AI12CS( 6) /   -.0000000263 1468847E0 /
      DATA AI12CS( 7) /   -.0000000038 3538039E0 /
      DATA AI12CS( 8) /   -.0000000005 5897433E0 /
      DATA AI12CS( 9) /   -.0000000000 1897495E0 /
      DATA AI12CS(10) /    .0000000000 3252602E0 /
      DATA AI12CS(11) /    .0000000000 1412580E0 /
      DATA AI12CS(12) /    .0000000000 0203564E0 /
      DATA AI12CS(13) /   -.0000000000 0071985E0 /
      DATA AI12CS(14) /   -.0000000000 0040836E0 /
      DATA AI12CS(15) /   -.0000000000 0002101E0 /
      DATA AI12CS(16) /    .0000000000 0004273E0 /
      DATA AI12CS(17) /    .0000000000 0001041E0 /
      DATA AI12CS(18) /   -.0000000000 0000382E0 /
      DATA AI12CS(19) /   -.0000000000 0000186E0 /
      DATA AI12CS(20) /    .0000000000 0000033E0 /
      DATA AI12CS(21) /    .0000000000 0000028E0 /
      DATA AI12CS(22) /   -.0000000000 0000003E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESI1E
      IF (FIRST) THEN
         NTI1 = INITS (BI1CS, 11, 0.1*R1MACH(3))
         NTAI1 = INITS (AI1CS, 21, 0.1*R1MACH(3))
         NTAI12 = INITS (AI12CS, 22, 0.1*R1MACH(3))
C
         XMIN = 2.0*R1MACH(1)
         XSML = SQRT (4.5*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0) GO TO 20
C
      BESI1E = 0.0
      IF (Y.EQ.0.0)  RETURN
C
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'BESI1E',
     +   'ABS(X) SO SMALL I1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) BESI1E = 0.5*X
      IF (Y.GT.XSML) BESI1E = X * (.875 + CSEVL(Y*Y/4.5-1., BI1CS,NTI1))
      BESI1E = EXP(-Y) * BESI1E
      RETURN
C
 20   IF (Y.LE.8.) BESI1E = (.375 + CSEVL ((48./Y-11.)/5., AI1CS, NTAI1)
     1  ) / SQRT(Y)
      IF (Y.GT.8.) BESI1E = (.375 + CSEVL (16./Y-1.0, AI12CS, NTAI12))
     1  / SQRT(Y)
      BESI1E = SIGN (BESI1E, X)
C
      RETURN
      END
