*DECK MPDIVI
      SUBROUTINE MPDIVI (X, IY, Z)
C***BEGIN PROLOGUE  MPDIVI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPDIVI-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Divides 'mp' X by the single-precision integer IY giving 'mp' Z.
C  This is much faster than division by an 'mp' number.
C
C  The arguments X(*) and Z(*), and the variable R in COMMON are all
C  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
C  for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPSTR, MPUNFL
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPDIVI
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), Z(*), RS, RE, R1, C, C2, B2
C***FIRST EXECUTABLE STATEMENT  MPDIVI
      RS = X(1)
      J = IY
      IF (J) 30, 10, 40
   10 WRITE (LUN, 20)
   20 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN CALL TO MPDIVI ***')
      GO TO 230
   30 J = -J
      RS = -RS
   40 RE = X(2)
C CHECK FOR ZERO DIVIDEND
      IF (RS.EQ.0) GO TO 120
C CHECK FOR DIVISION BY B
      IF (J.NE.B) GO TO 50
      CALL MPSTR (X, Z)
      IF (RE.LE.(-M)) GO TO 240
      Z(1) = RS
      Z(2) = RE - 1
      RETURN
C CHECK FOR DIVISION BY 1 OR -1
   50 IF (J.NE.1) GO TO 60
      CALL MPSTR (X, Z)
      Z(1) = RS
      RETURN
   60 C = 0
      I2 = T + 4
      I = 0
C IF J*B NOT REPRESENTABLE AS AN INTEGER HAVE TO SIMULATE
C LONG DIVISION.   ASSUME AT LEAST 16-BIT WORD.
      B2 = MAX(8*B,32767/B)
      IF (J.GE.B2) GO TO 130
C LOOK FOR FIRST NONZERO DIGIT IN QUOTIENT
   70 I = I + 1
      C = B*C
      IF (I.LE.T) C = C + X(I+2)
      R1 = C/J
      IF (R1) 210, 70, 80
C ADJUST EXPONENT AND GET T+4 DIGITS IN QUOTIENT
   80 RE = RE + 1 - I
      R(1) = R1
      C = B*(C - J*R1)
      KH = 2
      IF (I.GE.T) GO TO 100
      KH = 1 + T - I
      DO 90 K = 2, KH
      I = I + 1
      C = C + X(I+2)
      R(K) = C/J
   90 C = B*(C - J*R(K))
      IF (C.LT.0) GO TO 210
      KH = KH + 1
  100 DO 110 K = KH, I2
      R(K) = C/J
  110 C = B*(C - J*R(K))
      IF (C.LT.0) GO TO 210
C NORMALIZE AND ROUND RESULT
  120 CALL MPNZR (RS, RE, Z, 0)
      RETURN
C HERE NEED SIMULATED DOUBLE-PRECISION DIVISION
  130 C2 = 0
      J1 = J/B
      J2 = J - J1*B
      J11 = J1 + 1
C LOOK FOR FIRST NONZERO DIGIT
  140 I = I + 1
      C = B*C + C2
      C2 = 0
      IF (I.LE.T) C2 = X(I+2)
      IF (C-J1) 140, 150, 160
  150 IF (C2.LT.J2) GO TO 140
C COMPUTE T+4 QUOTIENT DIGITS
  160 RE = RE + 1 - I
      K = 1
      GO TO 180
C MAIN LOOP FOR LARGE ABS(IY) CASE
  170 K = K + 1
      IF (K.GT.I2) GO TO 120
      I = I + 1
C GET APPROXIMATE QUOTIENT FIRST
  180 IR = C/J11
C NOW REDUCE SO OVERFLOW DOES NOT OCCUR
      IQ = C - IR*J1
      IF (IQ.LT.B2) GO TO 190
C HERE IQ*B WOULD POSSIBLY OVERFLOW SO INCREASE IR
      IR = IR + 1
      IQ = IQ - J1
  190 IQ = IQ*B - IR*J2
      IF (IQ.GE.0) GO TO 200
C HERE IQ NEGATIVE SO IR WAS TOO LARGE
      IR = IR - 1
      IQ = IQ + J
  200 IF (I.LE.T) IQ = IQ + X(I+2)
      IQJ = IQ/J
C R(K) = QUOTIENT, C = REMAINDER
      R(K) = IQJ + IR
      C = IQ - J*IQJ
      IF (C.GE.0) GO TO 170
C CARRY NEGATIVE SO OVERFLOW MUST HAVE OCCURRED
  210 CALL MPCHK (1, 4)
      WRITE (LUN, 220)
  220 FORMAT (' *** INTEGER OVERFLOW IN MPDIVI, B TOO LARGE ***')
  230 CALL MPERR
      Z(1) = 0
      RETURN
C UNDERFLOW HERE
  240 CALL MPUNFL(Z)
      RETURN
      END
