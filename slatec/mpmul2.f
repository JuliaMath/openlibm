*DECK MPMUL2
      SUBROUTINE MPMUL2 (X, IY, Z, TRUNC)
C***BEGIN PROLOGUE  MPMUL2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPMUL2-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Multiplies 'mp' X by single-precision integer IY giving 'mp' Z.
C  Multiplication by 1 may be used to normalize a number even if some
C  digits are greater than B-1. Result is rounded if TRUNC.EQ.0,
C  otherwise truncated.
C
C  The arguments X(*) and Z(*), and the variable R in COMMON are all
C  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
C  for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPOVFL, MPSTR
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPMUL2
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), Z(*), TRUNC, RE, RS
      INTEGER C, C1, C2, RI, T1, T3, T4
C***FIRST EXECUTABLE STATEMENT  MPMUL2
      RS = X(1)
      IF (RS.EQ.0) GO TO 10
      J = IY
      IF (J) 20, 10, 50
C RESULT ZERO
   10 Z(1) = 0
      RETURN
   20 J = -J
      RS = -RS
C CHECK FOR MULTIPLICATION BY B
      IF (J.NE.B) GO TO 50
      IF (X(2).LT.M) GO TO 40
      CALL MPCHK (1, 4)
      WRITE (LUN, 30)
   30 FORMAT (' *** OVERFLOW OCCURRED IN MPMUL2 ***')
      CALL MPOVFL (Z)
      RETURN
   40 CALL MPSTR (X, Z)
      Z(1) = RS
      Z(2) = X(2) + 1
      RETURN
C SET EXPONENT TO EXPONENT(X) + 4
   50 RE = X(2) + 4
C FORM PRODUCT IN ACCUMULATOR
      C = 0
      T1 = T + 1
      T3 = T + 3
      T4 = T + 4
C IF J*B NOT REPRESENTABLE AS AN INTEGER WE HAVE TO SIMULATE
C DOUBLE-PRECISION MULTIPLICATION.
      IF (J.GE.MAX(8*B, 32767/B)) GO TO 110
      DO 60 IJ = 1, T
      I = T1 - IJ
      RI = J*X(I+2) + C
      C = RI/B
   60 R(I+4) = RI - B*C
C CHECK FOR INTEGER OVERFLOW
      IF (RI.LT.0) GO TO 130
C HAVE TO TREAT FIRST FOUR WORDS OF R SEPARATELY
      DO 70 IJ = 1, 4
      I = 5 - IJ
      RI = C
      C = RI/B
   70 R(I) = RI - B*C
      IF (C.EQ.0) GO TO 100
C HAVE TO SHIFT RIGHT HERE AS CARRY OFF END
   80 DO 90 IJ = 1, T3
      I = T4 - IJ
   90 R(I+1) = R(I)
      RI = C
      C = RI/B
      R(1) = RI - B*C
      RE = RE + 1
      IF (C) 130, 100, 80
C NORMALIZE AND ROUND OR TRUNCATE RESULT
  100 CALL MPNZR (RS, RE, Z, TRUNC)
      RETURN
C HERE J IS TOO LARGE FOR SINGLE-PRECISION MULTIPLICATION
  110 J1 = J/B
      J2 = J - J1*B
C FORM PRODUCT
      DO 120 IJ = 1, T4
      C1 = C/B
      C2 = C - B*C1
      I = T1 - IJ
      IX = 0
      IF (I.GT.0) IX = X(I+2)
      RI = J2*IX + C2
      IS = RI/B
      C = J1*IX + C1 + IS
  120 R(I+4) = RI - B*IS
      IF (C) 130, 100, 80
C CAN ONLY GET HERE IF INTEGER OVERFLOW OCCURRED
  130 CALL MPCHK (1, 4)
      WRITE (LUN, 140)
  140 FORMAT (' *** INTEGER OVERFLOW IN MPMUL2, B TOO LARGE ***')
      CALL MPERR
      GO TO 10
      END
