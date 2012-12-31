*DECK R9UPAK
      SUBROUTINE R9UPAK (X, Y, N)
C***BEGIN PROLOGUE  R9UPAK
C***PURPOSE  Unpack a floating point number X so that X = Y*2**N.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  A6B
C***TYPE      SINGLE PRECISION (R9UPAK-S, D9UPAK-D)
C***KEYWORDS  FNLIB, UNPACK
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   Unpack a floating point number X so that X = Y*2.0**N, where
C   0.5 .LE. ABS(Y) .LT. 1.0.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780701  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  R9UPAK
C***FIRST EXECUTABLE STATEMENT  R9UPAK
      ABSX = ABS(X)
      N = 0
      IF (X.EQ.0.0E0) GO TO 30
C
   10 IF (ABSX.GE.0.5E0) GO TO 20
      N = N-1
      ABSX = ABSX*2.0E0
      GO TO 10
C
   20 IF (ABSX.LT.1.0E0) GO TO 30
      N = N+1
      ABSX = ABSX*0.5E0
      GO TO 20
C
   30 Y = SIGN(ABSX,X)
      RETURN
C
      END
