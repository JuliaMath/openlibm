*DECK R9KNUS
      SUBROUTINE R9KNUS (XNU, X, BKNU, BKNU1, ISWTCH)
C***BEGIN PROLOGUE  R9KNUS
C***SUBSIDIARY
C***PURPOSE  Compute Bessel functions EXP(X)*K-SUB-XNU(X) and EXP(X)*
C            K-SUB-XNU+1(X) for 0.0 .LE. XNU .LT. 1.0.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B3
C***TYPE      SINGLE PRECISION (R9KNUS-S, D9KNUS-D)
C***KEYWORDS  BESSEL FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute Bessel functions EXP(X) * K-sub-XNU (X)  and
C EXP(X) * K-sub-XNU+1 (X) for 0.0 .LE. XNU .LT. 1.0 .
C
C Series for C0K        on the interval  0.          to  2.50000D-01
C                                        with weighted error   1.60E-17
C                                         log weighted error  16.79
C                               significant figures required  15.99
C                                    decimal places required  17.40
C
C Series for ZNU1       on the interval -7.00000D-01 to  0.
C                                        with weighted error   1.43E-17
C                                         log weighted error  16.85
C                               significant figures required  16.08
C                                    decimal places required  17.38
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, GAMMA, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  R9KNUS
      DIMENSION ALPHA(15), BETA(15), A(15), C0KCS(16), ZNU1CS(12)
      LOGICAL FIRST
      EXTERNAL GAMMA
      SAVE C0KCS, ZNU1CS, EULER, SQPI2, ALN2, NTC0K, NTZNU1,
     1 XNUSML, XSML, ALNSML, ALNBIG, ALNEPS, FIRST
      DATA C0KCS( 1) /    .0601830572 42626108E0 /
      DATA C0KCS( 2) /   -.1536487143 3017286E0 /
      DATA C0KCS( 3) /   -.0117511760 08210492E0 /
      DATA C0KCS( 4) /   -.0008524878 88919795E0 /
      DATA C0KCS( 5) /   -.0000613298 38767496E0 /
      DATA C0KCS( 6) /   -.0000044052 28124551E0 /
      DATA C0KCS( 7) /   -.0000003163 12467283E0 /
      DATA C0KCS( 8) /   -.0000000227 10719382E0 /
      DATA C0KCS( 9) /   -.0000000016 30564460E0 /
      DATA C0KCS(10) /   -.0000000001 17069392E0 /
      DATA C0KCS(11) /   -.0000000000 08405206E0 /
      DATA C0KCS(12) /   -.0000000000 00603466E0 /
      DATA C0KCS(13) /   -.0000000000 00043326E0 /
      DATA C0KCS(14) /   -.0000000000 00003110E0 /
      DATA C0KCS(15) /   -.0000000000 00000223E0 /
      DATA C0KCS(16) /   -.0000000000 00000016E0 /
      DATA ZNU1CS( 1) /    .2033067569 9419173E0 /
      DATA ZNU1CS( 2) /    .1400779334 1321977E0 /
      DATA ZNU1CS( 3) /    .0079167969 61001613E0 /
      DATA ZNU1CS( 4) /    .0003398011 82532104E0 /
      DATA ZNU1CS( 5) /    .0000117419 75688989E0 /
      DATA ZNU1CS( 6) /    .0000003393 57570612E0 /
      DATA ZNU1CS( 7) /    .0000000084 25941769E0 /
      DATA ZNU1CS( 8) /    .0000000001 83336677E0 /
      DATA ZNU1CS( 9) /    .0000000000 03549698E0 /
      DATA ZNU1CS(10) /    .0000000000 00061903E0 /
      DATA ZNU1CS(11) /    .0000000000 00000981E0 /
      DATA ZNU1CS(12) /    .0000000000 00000014E0 /
      DATA EULER / 0.5772156649 0153286E0 /
      DATA SQPI2 / 1.253314137 3155003E0 /
      DATA ALN2 / 0.693147180 55994531E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  R9KNUS
      IF (FIRST) THEN
         NTC0K = INITS (C0KCS, 16, 0.1*R1MACH(3))
         NTZNU1 = INITS (ZNU1CS, 12, 0.1*R1MACH(3))
C
         XNUSML = SQRT (R1MACH(3)/8.0)
         XSML = 0.1*R1MACH(3)
         ALNSML = LOG (R1MACH(1))
         ALNBIG = LOG (R1MACH(2))
         ALNEPS = LOG (0.1*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (XNU .LT. 0. .OR. XNU .GE. 1.0) CALL XERMSG ('SLATEC',
     +   'R9KNUS', 'XNU MUST BE GE 0 AND LT 1', 1, 2)
      IF (X .LE. 0.) CALL XERMSG ('SLATEC', 'R9KNUS', 'X MUST BE GT 0',
     +   2, 2)
C
      ISWTCH = 0
      IF (X.GT.2.0) GO TO 50
C
C X IS SMALL.  COMPUTE K-SUB-XNU (X) AND THE DERIVATIVE OF K-SUB-XNU (X)
C THEN FIND K-SUB-XNU+1 (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5)
C THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE
C ORDER (+NU).
C
      V = XNU
      IF (XNU.GT.0.5) V = 1.0 - XNU
C
C CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
      ALNZ = 2.0 * (LOG(X) - ALN2)
C
      IF (X.GT.XNU) GO TO 20
      IF (-0.5*XNU*ALNZ-ALN2-LOG(XNU) .GT. ALNBIG) CALL XERMSG
     +   ('SLATEC', 'R9KNUS', 'X SO SMALL BESSEL K-SUB-XNU OVERFLOWS',
     +   3, 2)
C
 20   VLNZ = V*ALNZ
      X2TOV = EXP (0.5*VLNZ)
      ZTOV = 0.0
      IF (VLNZ.GT.ALNSML) ZTOV = X2TOV**2
C
      A0 = 0.5*GAMMA(1.0+V)
      B0 = 0.5*GAMMA(1.0-V)
      C0 = -EULER
      IF (ZTOV.GT.0.5 .AND. V.GT.XNUSML) C0 = -0.75 +
     1  CSEVL ((8.0*V)*V-1., C0KCS, NTC0K)
C
      IF (ZTOV.LE.0.5) ALPHA(1) = (A0-ZTOV*B0)/V
      IF (ZTOV.GT.0.5) ALPHA(1) = C0 - ALNZ*(0.75 +
     1  CSEVL (VLNZ/0.35+1.0, ZNU1CS, NTZNU1))*B0
      BETA(1) = -0.5*(A0+ZTOV*B0)
C
      Z = 0.0
      IF (X.GT.XSML) Z = 0.25*X*X
      NTERMS = MAX (2.0, 11.0+(8.*ALNZ-25.19-ALNEPS)/(4.28-ALNZ))
      DO 30 I=2,NTERMS
        XI = I - 1
        A0 = A0/(XI*(XI-V))
        B0 = B0/(XI*(XI+V))
        ALPHA(I) = (ALPHA(I-1)+2.0*XI*A0)/(XI*(XI+V))
        BETA(I) = (XI-0.5*V)*ALPHA(I) - ZTOV*B0
 30   CONTINUE
C
      BKNU = ALPHA(NTERMS)
      BKNUD = BETA(NTERMS)
      DO 40 II=2,NTERMS
        I = NTERMS + 1 - II
        BKNU = ALPHA(I) + BKNU*Z
        BKNUD = BETA(I) + BKNUD*Z
 40   CONTINUE
C
      EXPX = EXP(X)
      BKNU = EXPX*BKNU/X2TOV
C
      IF (-0.5*(XNU+1.)*ALNZ-2.0*ALN2.GT.ALNBIG) ISWTCH = 1
      IF (ISWTCH.EQ.1) RETURN
      BKNUD = EXPX*BKNUD*2.0/(X2TOV*X)
C
      IF (XNU.LE.0.5) BKNU1 = V*BKNU/X - BKNUD
      IF (XNU.LE.0.5) RETURN
C
      BKNU0 = BKNU
      BKNU = -V*BKNU/X - BKNUD
      BKNU1 = 2.0*XNU*BKNU/X + BKNU0
      RETURN
C
C X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S
C RATIONAL EXPANSION.
C
 50   SQRTX = SQRT(X)
      IF (X.GT.1.0/XSML) GO TO 90
      AN = -1.56 + 4.0/X
      BN = -0.29 - 0.22/X
      NTERMS = MIN (15, MAX1 (3.0, AN+BN*ALNEPS))
C
      DO 80 INU=1,2
        XMU = 0.
        IF (INU.EQ.1 .AND. XNU.GT.XNUSML) XMU = (4.0*XNU)*XNU
        IF (INU.EQ.2) XMU = 4.0*(ABS(XNU)+1.)**2
C
        A(1) = 1.0 - XMU
        A(2) = 9.0 - XMU
        A(3) = 25.0 - XMU
        IF (A(2).EQ.0.) RESULT = SQPI2*(16.*X+XMU+7.)/(16.*X*SQRTX)
        IF (A(2).EQ.0.) GO TO 70
C
        ALPHA(1) = 1.0
        ALPHA(2) = (16.*X+A(2))/A(2)
        ALPHA(3) = ((768.*X+48.*A(3))*X + A(2)*A(3))/(A(2)*A(3))
C
        BETA(1) = 1.0
        BETA(2) = (16.*X+(XMU+7.))/A(2)
        BETA(3) = ((768.*X+48.*(XMU+23.))*X + ((XMU+62.)*XMU+129.))
     1    / (A(2)*A(3))
C
        IF (NTERMS.LT.4) GO TO 65
        DO 60 I=4,NTERMS
          N = I - 1
          X2N = 2*N - 1
C
          A(I) = (X2N+2.)**2 - XMU
          QQ = 16.*X2N/A(I)
          P1 = -X2N*(12*N*N-20*N-A(1))/((X2N-2.)*A(I)) - QQ*X
          P2 = (12*N*N-28*N+8-A(1))/A(I) - QQ*X
          P3 = -X2N*A(I-3)/((X2N-2.)*A(I))
C
          ALPHA(I) = -P1*ALPHA(I-1) - P2*ALPHA(I-2) - P3*ALPHA(I-3)
          BETA(I) = -P1*BETA(I-1) - P2*BETA(I-2) - P3*BETA(I-3)
 60     CONTINUE
C
 65     RESULT = SQPI2*BETA(NTERMS)/(SQRTX*ALPHA(NTERMS))
C
 70     IF (INU.EQ.1) BKNU = RESULT
        IF (INU.EQ.2) BKNU1 = RESULT
 80   CONTINUE
      RETURN
C
 90   BKNU = SQPI2/SQRTX
      BKNU1 = BKNU
      RETURN
C
      END
