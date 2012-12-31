*DECK R9LN2R
      FUNCTION R9LN2R (X)
C***BEGIN PROLOGUE  R9LN2R
C***SUBSIDIARY
C***PURPOSE  Evaluate LOG(1+X) from second order relative accuracy so
C            that LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      SINGLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so
C that    LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X)
C
C Series for LN21       on the interval -6.25000D-01 to  0.
C                                        with weighted error   2.49E-17
C                                         log weighted error  16.60
C                               significant figures required  15.87
C                                    decimal places required  17.31
C
C Series for LN22       on the interval  0.          to  8.12500D-01
C                                        with weighted error   1.42E-17
C                                         log weighted error  16.85
C                               significant figures required  15.95
C                                    decimal places required  17.50
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  R9LN2R
      REAL LN21CS(26), LN22CS(20)
      LOGICAL FIRST
      SAVE LN21CS, LN22CS, NTLN21, NTLN22, XMIN, XBIG, XMAX, FIRST
      DATA LN21CS( 1) /    .1811196251 3478810E0 /
      DATA LN21CS( 2) /   -.1562712319 2872463E0 /
      DATA LN21CS( 3) /    .0286763053 61557275E0 /
      DATA LN21CS( 4) /   -.0055586996 55948139E0 /
      DATA LN21CS( 5) /    .0011178976 65229983E0 /
      DATA LN21CS( 6) /   -.0002308050 89823279E0 /
      DATA LN21CS( 7) /    .0000485988 53341100E0 /
      DATA LN21CS( 8) /   -.0000103901 27388903E0 /
      DATA LN21CS( 9) /    .0000022484 56370739E0 /
      DATA LN21CS(10) /   -.0000004914 05927392E0 /
      DATA LN21CS(11) /    .0000001082 82565070E0 /
      DATA LN21CS(12) /   -.0000000240 25872763E0 /
      DATA LN21CS(13) /    .0000000053 62460047E0 /
      DATA LN21CS(14) /   -.0000000012 02995136E0 /
      DATA LN21CS(15) /    .0000000002 71078892E0 /
      DATA LN21CS(16) /   -.0000000000 61323562E0 /
      DATA LN21CS(17) /    .0000000000 13920858E0 /
      DATA LN21CS(18) /   -.0000000000 03169930E0 /
      DATA LN21CS(19) /    .0000000000 00723837E0 /
      DATA LN21CS(20) /   -.0000000000 00165700E0 /
      DATA LN21CS(21) /    .0000000000 00038018E0 /
      DATA LN21CS(22) /   -.0000000000 00008741E0 /
      DATA LN21CS(23) /    .0000000000 00002013E0 /
      DATA LN21CS(24) /   -.0000000000 00000464E0 /
      DATA LN21CS(25) /    .0000000000 00000107E0 /
      DATA LN21CS(26) /   -.0000000000 00000024E0 /
      DATA LN22CS( 1) /   -.2224253253 5020461E0 /
      DATA LN22CS( 2) /   -.0610471001 08078624E0 /
      DATA LN22CS( 3) /    .0074272350 09750394E0 /
      DATA LN22CS( 4) /   -.0009335018 26163697E0 /
      DATA LN22CS( 5) /    .0001200499 07687260E0 /
      DATA LN22CS( 6) /   -.0000157047 22952820E0 /
      DATA LN22CS( 7) /    .0000020818 74781051E0 /
      DATA LN22CS( 8) /   -.0000002789 19557764E0 /
      DATA LN22CS( 9) /    .0000000376 93558237E0 /
      DATA LN22CS(10) /   -.0000000051 30902896E0 /
      DATA LN22CS(11) /    .0000000007 02714117E0 /
      DATA LN22CS(12) /   -.0000000000 96748595E0 /
      DATA LN22CS(13) /    .0000000000 13381046E0 /
      DATA LN22CS(14) /   -.0000000000 01858102E0 /
      DATA LN22CS(15) /    .0000000000 00258929E0 /
      DATA LN22CS(16) /   -.0000000000 00036195E0 /
      DATA LN22CS(17) /    .0000000000 00005074E0 /
      DATA LN22CS(18) /   -.0000000000 00000713E0 /
      DATA LN22CS(19) /    .0000000000 00000100E0 /
      DATA LN22CS(20) /   -.0000000000 00000014E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  R9LN2R
      IF (FIRST) THEN
         EPS = R1MACH(3)
         NTLN21 = INITS (LN21CS, 26, 0.1*EPS)
         NTLN22 = INITS (LN22CS, 20, 0.1*EPS)
C
         XMIN = -1.0 + SQRT(R1MACH(4))
         SQEPS = SQRT(EPS)
         TXMAX = 6.0/SQEPS
         XMAX = TXMAX - (EPS*TXMAX**2 - 2.0*LOG(TXMAX)) /
     1                                              (2.0*EPS*TXMAX)
         TXBIG = 4.0/SQRT(SQEPS)
         XBIG = TXBIG - (SQEPS*TXBIG**2 - 2.0*LOG(TXBIG)) /
     1                                                (2.*SQEPS*TXBIG)
      ENDIF
      FIRST = .FALSE.
C
      IF (X.LT.(-0.625) .OR. X.GT.0.8125) GO TO 20
C
      IF (X.LT.0.0) R9LN2R = 0.375 + CSEVL (16.*X/5.+1.0, LN21CS,
     1  NTLN21)
      IF (X.GE.0.0) R9LN2R = 0.375 + CSEVL (32.*X/13.-1.0, LN22CS,
     1  NTLN22)
      RETURN
C
 20   IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'R9LN2R',
     +   'ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1', 1, 1)
      IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'R9LN2R',
     +   'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG', 3, 2)
      IF (X .GT. XBIG) CALL XERMSG ('SLATEC', 'R9LN2R',
     +   'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG', 2, 1)
C
      R9LN2R = (LOG(1.0+X) - X*(1.0-0.5*X) ) / X**3
      RETURN
C
      END
