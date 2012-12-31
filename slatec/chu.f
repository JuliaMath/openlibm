*DECK CHU
      FUNCTION CHU (A, B, X)
C***BEGIN PROLOGUE  CHU
C***PURPOSE  Compute the logarithmic confluent hypergeometric function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C11
C***TYPE      SINGLE PRECISION (CHU-S, DCHU-D)
C***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CHU computes the logarithmic confluent hypergeometric function,
C U(A,B,X).
C
C Input Parameters:
C       A   real
C       B   real
C       X   real and positive
C
C This routine is not valid when 1+A-B is close to zero if X is small.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  EXPREL, GAMMA, GAMR, POCH, POCH1, R1MACH, R9CHU,
C                    XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  CHU
      EXTERNAL GAMMA
      SAVE PI, EPS
      DATA PI / 3.1415926535 8979324 E0 /
      DATA EPS / 0.0 /
C***FIRST EXECUTABLE STATEMENT  CHU
      IF (EPS.EQ.0.0) EPS = R1MACH(3)
C
      IF (X .EQ. 0.0) CALL XERMSG ('SLATEC', 'CHU',
     +   'X IS ZERO SO CHU IS INFINITE', 1, 2)
      IF (X .LT. 0.0) CALL XERMSG ('SLATEC', 'CHU',
     +   'X IS NEGATIVE, USE CCHU', 2, 2)
C
      IF (MAX(ABS(A),1.0)*MAX(ABS(1.0+A-B),1.0).LT.0.99*ABS(X))
     1  GO TO 120
C
C THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL
C APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE.
C
      IF (ABS(1.0+A-B) .LT. SQRT(EPS)) CALL XERMSG ('SLATEC', 'CHU',
     +   'ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO FOR SMALL X', 10, 2)
C
      AINTB = AINT(B+0.5)
      IF (B.LT.0.0) AINTB = AINT(B-0.5)
      BEPS = B - AINTB
      N = AINTB
C
      ALNX = LOG(X)
      XTOEPS = EXP(-BEPS*ALNX)
C
C EVALUATE THE FINITE SUM.     -----------------------------------------
C
      IF (N.GE.1) GO TO 40
C
C CONSIDER THE CASE B .LT. 1.0 FIRST.
C
      SUM = 1.0
      IF (N.EQ.0) GO TO 30
C
      T = 1.0
      M = -N
      DO 20 I=1,M
        XI1 = I - 1
        T = T*(A+XI1)*X/((B+XI1)*(XI1+1.0))
        SUM = SUM + T
 20   CONTINUE
C
 30   SUM = POCH(1.0+A-B, -A) * SUM
      GO TO 70
C
C NOW CONSIDER THE CASE B .GE. 1.0.
C
 40   SUM = 0.0
      M = N - 2
      IF (M.LT.0) GO TO 70
      T = 1.0
      SUM = 1.0
      IF (M.EQ.0) GO TO 60
C
      DO 50 I=1,M
        XI = I
        T = T * (A-B+XI)*X/((1.0-B+XI)*XI)
        SUM = SUM + T
 50   CONTINUE
C
 60   SUM = GAMMA(B-1.0) * GAMR(A) * X**(1-N) * XTOEPS * SUM
C
C NOW EVALUATE THE INFINITE SUM.     -----------------------------------
C
 70   ISTRT = 0
      IF (N.LT.1) ISTRT = 1 - N
      XI = ISTRT
C
      FACTOR = (-1.0)**N * GAMR(1.0+A-B) * X**ISTRT
      IF (BEPS.NE.0.0) FACTOR = FACTOR * BEPS*PI/SIN(BEPS*PI)
C
      POCHAI = POCH (A, XI)
      GAMRI1 = GAMR (XI+1.0)
      GAMRNI = GAMR (AINTB+XI)
      B0 = FACTOR * POCH(A,XI-BEPS) * GAMRNI * GAMR(XI+1.0-BEPS)
C
      IF (ABS(XTOEPS-1.0).GT.0.5) GO TO 90
C
C X**(-BEPS) IS CLOSE TO 1.0, SO WE MUST BE CAREFUL IN EVALUATING
C THE DIFFERENCES
C
      PCH1AI = POCH1 (A+XI, -BEPS)
      PCH1I = POCH1 (XI+1.0-BEPS, BEPS)
      C0 = FACTOR * POCHAI * GAMRNI * GAMRI1 * (
     1  -POCH1(B+XI, -BEPS) + PCH1AI - PCH1I + BEPS*PCH1AI*PCH1I )
C
C XEPS1 = (1.0 - X**(-BEPS)) / BEPS
      XEPS1 = ALNX * EXPREL(-BEPS*ALNX)
C
      CHU = SUM + C0 + XEPS1*B0
      XN = N
      DO 80 I=1,1000
        XI = ISTRT + I
        XI1 = ISTRT + I - 1
        B0 = (A+XI1-BEPS)*B0*X/((XN+XI1)*(XI-BEPS))
        C0 = (A+XI1)*C0*X/((B+XI1)*XI) - ((A-1.0)*(XN+2.*XI-1.0)
     1    + XI*(XI-BEPS)) * B0/(XI*(B+XI1)*(A+XI1-BEPS))
        T = C0 + XEPS1*B0
        CHU = CHU + T
        IF (ABS(T).LT.EPS*ABS(CHU)) GO TO 130
 80   CONTINUE
      CALL XERMSG ('SLATEC', 'CHU',
     +   'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
C
C X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD
C FORMULATION IS STABLE.
C
 90   A0 = FACTOR * POCHAI * GAMR(B+XI) * GAMRI1 / BEPS
      B0 = XTOEPS*B0/BEPS
C
      CHU = SUM + A0 - B0
      DO 100 I=1,1000
        XI = ISTRT + I
        XI1 = ISTRT + I - 1
        A0 = (A+XI1)*A0*X/((B+XI1)*XI)
        B0 = (A+XI1-BEPS)*B0*X/((AINTB+XI1)*(XI-BEPS))
        T = A0 - B0
        CHU = CHU + T
        IF (ABS(T).LT.EPS*ABS(CHU)) GO TO 130
 100  CONTINUE
      CALL XERMSG ('SLATEC', 'CHU',
     +   'NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING SERIES', 3, 2)
C
C USE LUKE-S RATIONAL APPROX IN THE ASYMPTOTIC REGION.
C
 120  CHU = X**(-A) * R9CHU(A, B, X)
C
 130  RETURN
      END
