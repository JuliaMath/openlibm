*DECK BIE
      FUNCTION BIE (X)
C***BEGIN PROLOGUE  BIE
C***PURPOSE  Calculate the Bairy function for a negative argument and an
C            exponentially scaled Bairy function for a non-negative
C            argument.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10D
C***TYPE      SINGLE PRECISION (BIE-S, DBIE-D)
C***KEYWORDS  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate BI(X) for X .LE. 0  and  BI(X)*EXP(ZETA)  where
C ZETA = 2/3 * X**(3/2)  for X .GE. 0.0
C
C Series for BIF        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.88E-19
C                                         log weighted error  18.72
C                               significant figures required  17.74
C                                    decimal places required  19.20
C
C Series for BIG        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   2.61E-17
C                                         log weighted error  16.58
C                               significant figures required  15.17
C                                    decimal places required  17.03
C
C Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00
C                                        with weighted error   1.11E-17
C                                         log weighted error  16.95
C                        approx significant figures required  16.5
C                                    decimal places required  17.45
C
C Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00
C                                        with weighted error   1.19E-18
C                                         log weighted error  17.92
C                        approx significant figures required  17.2
C                                    decimal places required  18.42
C
C Series for BIP        on the interval  1.25000D-01 to  3.53553D-01
C                                        with weighted error   1.91E-17
C                                         log weighted error  16.72
C                               significant figures required  15.35
C                                    decimal places required  17.41
C
C Series for BIP2       on the interval  0.          to  1.25000D-01
C                                        with weighted error   1.05E-18
C                                         log weighted error  17.98
C                               significant figures required  16.74
C                                    decimal places required  18.71
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  BIE
      LOGICAL FIRST
      DIMENSION BIFCS(9), BIGCS(8), BIF2CS(10), BIG2CS(10), BIPCS(24),
     1  BIP2CS(29)
      SAVE BIFCS, BIGCS, BIF2CS, BIG2CS, BIPCS, BIP2CS, ATR, BTR,
     1 NBIF, NBIG, NBIF2, NBIG2, NBIP, NBIP2, X3SML, X32SML, XBIG, FIRST
      DATA BIFCS( 1) /   -.0167302164 7198664948E0 /
      DATA BIFCS( 2) /    .1025233583 424944561E0 /
      DATA BIFCS( 3) /    .0017083092 5073815165E0 /
      DATA BIFCS( 4) /    .0000118625 4546774468E0 /
      DATA BIFCS( 5) /    .0000000449 3290701779E0 /
      DATA BIFCS( 6) /    .0000000001 0698207143E0 /
      DATA BIFCS( 7) /    .0000000000 0017480643E0 /
      DATA BIFCS( 8) /    .0000000000 0000020810E0 /
      DATA BIFCS( 9) /    .0000000000 0000000018E0 /
      DATA BIGCS( 1) /    .0224662232 4857452E0 /
      DATA BIGCS( 2) /    .0373647754 5301955E0 /
      DATA BIGCS( 3) /    .0004447621 8957212E0 /
      DATA BIGCS( 4) /    .0000024708 0756363E0 /
      DATA BIGCS( 5) /    .0000000079 1913533E0 /
      DATA BIGCS( 6) /    .0000000000 1649807E0 /
      DATA BIGCS( 7) /    .0000000000 0002411E0 /
      DATA BIGCS( 8) /    .0000000000 0000002E0 /
      DATA BIF2CS( 1) /   0.0998457269 3816041E0 /
      DATA BIF2CS( 2) /    .4786249778 63005538E0 /
      DATA BIF2CS( 3) /    .0251552119 604330118E0 /
      DATA BIF2CS( 4) /    .0005820693 885232645E0 /
      DATA BIF2CS( 5) /    .0000074997 659644377E0 /
      DATA BIF2CS( 6) /    .0000000613 460287034E0 /
      DATA BIF2CS( 7) /    .0000000003 462753885E0 /
      DATA BIF2CS( 8) /    .0000000000 014288910E0 /
      DATA BIF2CS( 9) /    .0000000000 000044962E0 /
      DATA BIF2CS(10) /    .0000000000 000000111E0 /
      DATA BIG2CS( 1) /    .0333056621 45514340E0 /
      DATA BIG2CS( 2) /    .1613092151 23197068E0 /
      DATA BIG2CS( 3) /    .0063190073 096134286E0 /
      DATA BIG2CS( 4) /    .0001187904 568162517E0 /
      DATA BIG2CS( 5) /    .0000013045 345886200E0 /
      DATA BIG2CS( 6) /    .0000000093 741259955E0 /
      DATA BIG2CS( 7) /    .0000000000 474580188E0 /
      DATA BIG2CS( 8) /    .0000000000 001783107E0 /
      DATA BIG2CS( 9) /    .0000000000 000005167E0 /
      DATA BIG2CS(10) /    .0000000000 000000011E0 /
      DATA BIPCS( 1) /   -.0832204747 7943447E0 /
      DATA BIPCS( 2) /    .0114611892 7371174E0 /
      DATA BIPCS( 3) /    .0004289644 0718911E0 /
      DATA BIPCS( 4) /   -.0001490663 9379950E0 /
      DATA BIPCS( 5) /   -.0000130765 9726787E0 /
      DATA BIPCS( 6) /    .0000063275 9839610E0 /
      DATA BIPCS( 7) /   -.0000004222 6696982E0 /
      DATA BIPCS( 8) /   -.0000001914 7186298E0 /
      DATA BIPCS( 9) /    .0000000645 3106284E0 /
      DATA BIPCS(10) /   -.0000000078 4485467E0 /
      DATA BIPCS(11) /   -.0000000009 6077216E0 /
      DATA BIPCS(12) /    .0000000007 0004713E0 /
      DATA BIPCS(13) /   -.0000000001 7731789E0 /
      DATA BIPCS(14) /    .0000000000 2272089E0 /
      DATA BIPCS(15) /    .0000000000 0165404E0 /
      DATA BIPCS(16) /   -.0000000000 0185171E0 /
      DATA BIPCS(17) /    .0000000000 0059576E0 /
      DATA BIPCS(18) /   -.0000000000 0012194E0 /
      DATA BIPCS(19) /    .0000000000 0001334E0 /
      DATA BIPCS(20) /    .0000000000 0000172E0 /
      DATA BIPCS(21) /   -.0000000000 0000145E0 /
      DATA BIPCS(22) /    .0000000000 0000049E0 /
      DATA BIPCS(23) /   -.0000000000 0000011E0 /
      DATA BIPCS(24) /    .0000000000 0000001E0 /
      DATA BIP2CS( 1) /   -.1135967375 85988679E0 /
      DATA BIP2CS( 2) /    .0041381473 947881595E0 /
      DATA BIP2CS( 3) /    .0001353470 622119332E0 /
      DATA BIP2CS( 4) /    .0000104273 166530153E0 /
      DATA BIP2CS( 5) /    .0000013474 954767849E0 /
      DATA BIP2CS( 6) /    .0000001696 537405438E0 /
      DATA BIP2CS( 7) /   -.0000000100 965008656E0 /
      DATA BIP2CS( 8) /   -.0000000167 291194937E0 /
      DATA BIP2CS( 9) /   -.0000000045 815364485E0 /
      DATA BIP2CS(10) /    .0000000003 736681366E0 /
      DATA BIP2CS(11) /    .0000000005 766930320E0 /
      DATA BIP2CS(12) /    .0000000000 621812650E0 /
      DATA BIP2CS(13) /   -.0000000000 632941202E0 /
      DATA BIP2CS(14) /   -.0000000000 149150479E0 /
      DATA BIP2CS(15) /    .0000000000 078896213E0 /
      DATA BIP2CS(16) /    .0000000000 024960513E0 /
      DATA BIP2CS(17) /   -.0000000000 012130075E0 /
      DATA BIP2CS(18) /   -.0000000000 003740493E0 /
      DATA BIP2CS(19) /    .0000000000 002237727E0 /
      DATA BIP2CS(20) /    .0000000000 000474902E0 /
      DATA BIP2CS(21) /   -.0000000000 000452616E0 /
      DATA BIP2CS(22) /   -.0000000000 000030172E0 /
      DATA BIP2CS(23) /    .0000000000 000091058E0 /
      DATA BIP2CS(24) /   -.0000000000 000009814E0 /
      DATA BIP2CS(25) /   -.0000000000 000016429E0 /
      DATA BIP2CS(26) /    .0000000000 000005533E0 /
      DATA BIP2CS(27) /    .0000000000 000002175E0 /
      DATA BIP2CS(28) /   -.0000000000 000001737E0 /
      DATA BIP2CS(29) /   -.0000000000 000000010E0 /
      DATA ATR / 8.750690570 8484345 E0 /
      DATA BTR / -2.093836321 356054 E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BIE
      IF (FIRST) THEN
         ETA = 0.1*R1MACH(3)
         NBIF = INITS (BIFCS, 9, ETA)
         NBIG = INITS (BIGCS, 8, ETA)
         NBIF2 = INITS (BIF2CS, 10, ETA)
         NBIG2 = INITS (BIG2CS, 10, ETA)
         NBIP  = INITS (BIPCS , 24, ETA)
         NBIP2 = INITS (BIP2CS, 29, ETA)
C
         X3SML = ETA**0.3333
         X32SML = 1.3104*X3SML**2
         XBIG = R1MACH(2)**0.6666
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GE.(-1.0)) GO TO 20
      CALL R9AIMP (X, XM, THETA)
      BIE = XM * SIN(THETA)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 30
      Z = 0.0
      IF (ABS(X).GT.X3SML) Z = X**3
      BIE = 0.625 + CSEVL (Z, BIFCS, NBIF) + X*(0.4375 +
     1  CSEVL (Z, BIGCS, NBIG))
      IF (X.GT.X32SML) BIE = BIE * EXP(-2.0*X*SQRT(X)/3.0)
      RETURN
C
 30   IF (X.GT.2.0) GO TO 40
      Z = (2.0*X**3 - 9.0) / 7.0
      BIE = EXP(-2.0*X*SQRT(X)/3.0) * (1.125 + CSEVL (Z, BIF2CS, NBIF2)
     1  + X*(0.625 + CSEVL (Z, BIG2CS, NBIG2)) )
      RETURN
C
 40   IF (X.GT.4.0) GO TO 50
      SQRTX = SQRT(X)
      Z = ATR/(X*SQRTX) + BTR
      BIE = (0.625 + CSEVL (Z, BIPCS, NBIP)) / SQRT(SQRTX)
      RETURN
C
 50   SQRTX = SQRT(X)
      Z = -1.0
      IF (X.LT.XBIG) Z = 16.0/(X*SQRTX) - 1.0
      BIE = (0.625 + CSEVL (Z, BIP2CS, NBIP2))/SQRT(SQRTX)
      RETURN
C
      END
