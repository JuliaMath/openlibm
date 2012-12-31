*DECK MPADD3
      SUBROUTINE MPADD3 (X, Y, S, MED, RE)
C***BEGIN PROLOGUE  MPADD3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPADD3-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   Called by MPADD2; does inner loops of addition
C
C   The arguments X(*) and Y(*) and the variable R in COMMON are all
C   INTEGER arrays of size 30.  See the comments in the routine MPBLAS
C   for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPADD3
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), Y(*), S, RE, C, TED
C***FIRST EXECUTABLE STATEMENT  MPADD3
      TED = T + MED
      I2 = T + 4
      I = I2
      C = 0
C CLEAR GUARD DIGITS TO RIGHT OF X DIGITS
   10 IF (I.LE.TED) GO TO 20
      R(I) = 0
      I = I - 1
      GO TO 10
   20 IF (S.LT.0) GO TO 130
C HERE DO ADDITION, EXPONENT(Y) .GE. EXPONENT(X)
      IF (I.LT.T) GO TO 40
   30 J = I - MED
      R(I) = X(J+2)
      I = I - 1
      IF (I.GT.T) GO TO 30
   40 IF (I.LE.MED) GO TO 60
      J = I - MED
      C = Y(I+2) + X(J+2) + C
      IF (C.LT.B) GO TO 50
C CARRY GENERATED HERE
      R(I) = C - B
      C = 1
      I = I - 1
      GO TO 40
C NO CARRY GENERATED HERE
   50 R(I) = C
      C = 0
      I = I - 1
      GO TO 40
   60 IF (I.LE.0) GO TO 90
      C = Y(I+2) + C
      IF (C.LT.B) GO TO 70
      R(I) = 0
      C = 1
      I = I - 1
      GO TO 60
   70 R(I) = C
      I = I - 1
C NO CARRY POSSIBLE HERE
   80 IF (I.LE.0) RETURN
      R(I) = Y(I+2)
      I = I - 1
      GO TO 80
   90 IF (C.EQ.0) RETURN
C MUST SHIFT RIGHT HERE AS CARRY OFF END
      I2P = I2 + 1
      DO 100 J = 2, I2
      I = I2P - J
  100 R(I+1) = R(I)
      R(1) = 1
      RE = RE + 1
      RETURN
C HERE DO SUBTRACTION, ABS(Y) .GT. ABS(X)
  110 J = I - MED
      R(I) = C - X(J+2)
      C = 0
      IF (R(I).GE.0) GO TO 120
C BORROW GENERATED HERE
      C = -1
      R(I) = R(I) + B
  120 I = I - 1
  130 IF (I.GT.T) GO TO 110
  140 IF (I.LE.MED) GO TO 160
      J = I - MED
      C = Y(I+2) + C - X(J+2)
      IF (C.GE.0) GO TO 150
C BORROW GENERATED HERE
      R(I) = C + B
      C = -1
      I = I - 1
      GO TO 140
C NO BORROW GENERATED HERE
  150 R(I) = C
      C = 0
      I = I - 1
      GO TO 140
  160 IF (I.LE.0) RETURN
      C = Y(I+2) + C
      IF (C.GE.0) GO TO 70
      R(I) = C + B
      C = -1
      I = I - 1
      GO TO 160
      END
