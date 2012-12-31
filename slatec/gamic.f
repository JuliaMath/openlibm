*DECK GAMIC
      REAL FUNCTION GAMIC (A, X)
C***BEGIN PROLOGUE  GAMIC
C***PURPOSE  Calculate the complementary incomplete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      SINGLE PRECISION (GAMIC-S, DGAMIC-D)
C***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   Evaluate the complementary incomplete gamma function
C
C   GAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  .
C
C   GAMIC is evaluated for arbitrary real values of A and for non-
C   negative values of X (even though GAMIC is defined for X .LT.
C   0.0), except that for X = 0 and A .LE. 0.0, GAMIC is undefined.
C
C   GAMIC, A, and X are REAL.
C
C   A slight deterioration of 2 or 3 digits accuracy will occur when
C   GAMIC is very large or very small in absolute value, because log-
C   arithmic variables are used.  Also, if the parameter A is very close
C   to a negative integer (but not a negative integer), there is a loss
C   of accuracy, which is reported if the result is less than half
C   machine precision.
C
C***REFERENCES  W. Gautschi, A computational procedure for incomplete
C                 gamma functions, ACM Transactions on Mathematical
C                 Software 5, 4 (December 1979), pp. 466-481.
C               W. Gautschi, Incomplete gamma functions, Algorithm 542,
C                 ACM Transactions on Mathematical Software 5, 4
C                 (December 1979), pp. 482-489.
C***ROUTINES CALLED  ALGAMS, ALNGAM, R1MACH, R9GMIC, R9GMIT, R9LGIC,
C                    R9LGIT, XERCLR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  GAMIC
      LOGICAL FIRST
      SAVE EPS, SQEPS, ALNEPS, BOT, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  GAMIC
      IF (FIRST) THEN
         EPS = 0.5*R1MACH(3)
         SQEPS = SQRT(R1MACH(4))
         ALNEPS = -LOG(R1MACH(3))
         BOT = LOG(R1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 0.0) CALL XERMSG ('SLATEC', 'GAMIC', 'X IS NEGATIVE',
     +   2, 2)
C
      IF (X.GT.0.0) GO TO 20
      IF (A .LE. 0.0) CALL XERMSG ('SLATEC', 'GAMIC',
     +   'X = 0 AND A LE 0 SO GAMIC IS UNDEFINED', 3, 2)
C
      GAMIC = EXP (ALNGAM(A+1.0) - LOG(A))
      RETURN
C
 20   ALX = LOG(X)
      SGA = 1.0
      IF (A.NE.0.0) SGA = SIGN (1.0, A)
      MA = A + 0.5*SGA
      AEPS = A - MA
C
      IZERO = 0
      IF (X.GE.1.0) GO TO 60
C
      IF (A.GT.0.5 .OR. ABS(AEPS).GT.0.001) GO TO 50
      FM = -MA
      E = 2.0
      IF (FM.GT.1.0) E = 2.0*(FM+2.0)/(FM*FM-1.0)
      E = E - ALX*X**(-0.001)
      IF (E*ABS(AEPS).GT.EPS) GO TO 50
C
      GAMIC = R9GMIC (A, X, ALX)
      RETURN
C
 50   CALL ALGAMS (A+1.0, ALGAP1, SGNGAM)
      GSTAR = R9GMIT (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0.0) IZERO = 1
      IF (GSTAR.NE.0.0) ALNGS = LOG (ABS(GSTAR))
      IF (GSTAR.NE.0.0) SGNGS = SIGN (1.0, GSTAR)
      GO TO 70
C
 60   IF (A.LT.X) GAMIC = EXP (R9LGIC(A, X, ALX))
      IF (A.LT.X) RETURN
C
      SGNGAM = 1.0
      ALGAP1 = ALNGAM (A+1.0)
      SGNGS = 1.0
      ALNGS = R9LGIT (A, X, ALGAP1)
C
C EVALUATION OF GAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
C
 70   H = 1.0
      IF (IZERO.EQ.1) GO TO 80
C
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 90
      IF (T.GT.(-ALNEPS)) H = 1.0 - SGNGS*EXP(T)
C
      IF (ABS(H).LT.SQEPS) CALL XERCLR
      IF (ABS(H) .LT. SQEPS) CALL XERMSG ('SLATEC', 'GAMIC',
     +   'RESULT LT HALF PRECISION', 1, 1)
C
 80   SGNG = SIGN (1.0, H) * SGA * SGNGAM
      T = LOG(ABS(H)) + ALGAP1 - LOG(ABS(A))
      IF (T.LT.BOT) CALL XERCLR
      GAMIC = SGNG * EXP(T)
      RETURN
C
 90   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - LOG(ABS(A))
      IF (T.LT.BOT) CALL XERCLR
      GAMIC = SGNG * EXP(T)
      RETURN
C
      END
