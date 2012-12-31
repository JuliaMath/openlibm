*DECK MPNZR
      SUBROUTINE MPNZR (RS, RE, Z, TRUNC)
C***BEGIN PROLOGUE  MPNZR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPNZR-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Modified for use with BLAS.  Blank COMMON changed to named COMMON.
C  Assumes long (i.e. (t+4)-DIGIT) fraction in R, sign = RS, exponent
C  = RE.  Normalizes, and returns 'mp' result in Z. Integer arguments
C  RS and RE are not preserved. R*-rounding is used if TRUNC.EQ.0
C
C  The argument Z(*) and the variable R in COMMON are INTEGER arrays
C  of size 30.  See the comments in the routine MPBLAS for the reason
C  for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPERR, MPOVFL, MPUNFL
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPNZR
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, Z(*), RE, RS, TRUNC, B2
C***FIRST EXECUTABLE STATEMENT  MPNZR
      I2 = T + 4
      IF (RS.NE.0) GO TO 20
C STORE ZERO IN Z
   10 Z(1) = 0
      RETURN
C CHECK THAT SIGN = +-1
   20 IF (ABS(RS).LE.1) GO TO 40
      WRITE (LUN, 30)
   30 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN CALL TO MPNZR,',
     1        ' POSSIBLE OVERWRITING PROBLEM ***')
      CALL MPERR
      GO TO 10
C LOOK FOR FIRST NONZERO DIGIT
   40 DO 50 I = 1, I2
      IS = I - 1
      IF (R(I).GT.0) GO TO 60
   50 CONTINUE
C FRACTION ZERO
      GO TO 10
   60 IF (IS.EQ.0) GO TO 90
C NORMALIZE
      RE = RE - IS
      I2M = I2 - IS
      DO 70 J = 1, I2M
      K = J + IS
   70 R(J) = R(K)
      I2P = I2M + 1
      DO 80 J = I2P, I2
   80 R(J) = 0
C CHECK TO SEE IF TRUNCATION IS DESIRED
   90 IF (TRUNC.NE.0) GO TO 150
C SEE IF ROUNDING NECESSARY
C TREAT EVEN AND ODD BASES DIFFERENTLY
      B2 = B/2
      IF ((2*B2).NE.B) GO TO 130
C B EVEN.  ROUND IF R(T+1).GE.B2 UNLESS R(T) ODD AND ALL ZEROS
C AFTER R(T+2).
      IF (R(T+1) - B2) 150, 100, 110
  100 IF (MOD(R(T),2).EQ.0) GO TO 110
      IF ((R(T+2)+R(T+3)+R(T+4)).EQ.0) GO TO 150
C ROUND
  110 DO 120 J = 1, T
      I = T + 1 - J
      R(I) = R(I) + 1
      IF (R(I).LT.B) GO TO 150
  120 R(I) = 0
C EXCEPTIONAL CASE, ROUNDED UP TO .10000...
      RE = RE + 1
      R(1) = 1
      GO TO 150
C ODD BASE, ROUND IF R(T+1)... .GT. 1/2
  130 DO 140 I = 1, 4
      IT = T + I
      IF (R(IT) - B2) 150, 140, 110
  140 CONTINUE
C CHECK FOR OVERFLOW
  150 IF (RE.LE.M) GO TO 170
      WRITE (LUN, 160)
  160 FORMAT (' *** OVERFLOW OCCURRED IN MPNZR ***')
      CALL MPOVFL (Z)
      RETURN
C CHECK FOR UNDERFLOW
  170 IF (RE.LT.(-M)) GO TO 190
C STORE RESULT IN Z
      Z(1) = RS
      Z(2) = RE
      DO 180 I = 1, T
  180 Z(I+2) = R(I)
      RETURN
C UNDERFLOW HERE
  190 CALL MPUNFL (Z)
      RETURN
      END
