*DECK POCH
      FUNCTION POCH (A, X)
C***BEGIN PROLOGUE  POCH
C***PURPOSE  Evaluate a generalization of Pochhammer's symbol.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C1, C7A
C***TYPE      SINGLE PRECISION (POCH-S, DPOCH-D)
C***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate a generalization of Pochhammer's symbol
C (A)-sub-X = GAMMA(A+X)/GAMMA(A).  For X a non-negative integer,
C POCH(A,X) is just Pochhammer's symbol.  A and X are single precision.
C This is a preliminary version.  Error handling when POCH(A,X) is
C less than half precision is probably incorrect.  Grossly incorrect
C arguments are not handled properly.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALGAMS, ALNREL, FAC, GAMMA, GAMR, R9LGMC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  POCH
      EXTERNAL GAMMA
      SAVE PI
      DATA PI / 3.1415926535 89793238 E0 /
C***FIRST EXECUTABLE STATEMENT  POCH
      AX = A + X
      IF (AX.GT.0.0) GO TO 30
      IF (AINT(AX).NE.AX) GO TO 30
C
      IF (A .GT. 0.0 .OR. AINT(A) .NE. A) CALL XERMSG ('SLATEC', 'POCH',
     +   'A+X IS NON-POSITIVE INTEGER BUT A IS NOT', 2, 2)
C
C WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS.
C
      POCH = 1.0
      IF (X.EQ.0.0) RETURN
C
      N = X
      IF (MIN(A+X,A).LT.(-20.0)) GO TO 20
C
      POCH = (-1.0)**N * FAC(-INT(A))/FAC(-INT(A)-N)
      RETURN
C
 20   POCH = (-1.0)**N * EXP ((A-0.5)*ALNREL(X/(A-1.0))
     1  + X*LOG(-A+1.0-X) - X + R9LGMC(-A+1.) - R9LGMC(-A-X+1.) )
      RETURN
C
C HERE WE KNOW A+X IS NOT ZERO OR A NEGATIVE INTEGER.
C
 30   POCH = 0.0
      IF (A.LE.0.0 .AND. AINT(A).EQ.A) RETURN
C
      N = ABS(X)
      IF (REAL(N).NE.X .OR. N.GT.20) GO TO 50
C
C X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE.
C
      POCH = 1.0
      IF (N.EQ.0) RETURN
      DO 40 I=1,N
        POCH = POCH * (A+I-1)
 40   CONTINUE
      RETURN
C
 50   ABSAX = ABS(A+X)
      ABSA = ABS(A)
      IF (MAX(ABSAX,ABSA).GT.20.0) GO TO 60
      POCH = GAMMA(A+X)*GAMR(A)
      RETURN
C
 60   IF (ABS(X).GT.0.5*ABSA) GO TO 70
C
C HERE ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS,
C A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE
C GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) *
C SIN(PI*A)/SIN(PI*(A+X))
C
      B = A
      IF (B.LT.0.0) B = -A - X + 1.0
      POCH = EXP ((B-0.5)*ALNREL(X/B) + X*LOG(B+X) - X +
     1  R9LGMC(B+X) - R9LGMC(B) )
      IF (A.LT.0.0 .AND. POCH.NE.0.0) POCH = POCH/(COS(PI*X) +
     1  COT(PI*A)*SIN(PI*X))
      RETURN
C
 70   CALL ALGAMS (A+X, ALNGAX, SGNGAX)
      CALL ALGAMS (A, ALNGA, SGNGA)
      POCH = SGNGAX * SGNGA * EXP(ALNGAX-ALNGA)
C
      RETURN
      END
