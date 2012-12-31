*DECK MPMUL
      SUBROUTINE MPMUL (X, Y, Z)
C***BEGIN PROLOGUE  MPMUL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPMUL-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Multiplies X and Y, returning result in Z, for 'mp' X, Y and Z.
C  The simple o(t**2) algorithm is used, with four guard digits and
C  R*-rounding. Advantage is taken of zero digits in X, but not in Y.
C  Asymptotically faster algorithms are known (see Knuth, VOL. 2),
C  but are difficult to implement in FORTRAN in an efficient and
C  machine-independent manner. In comments to other 'mp' routines,
C  M(t) is the time to perform t-digit 'mp' multiplication. Thus
C  M(t) = o(t**2) with the present version of MPMUL, but
C  M(t) = o(t.log(t).log(log(t))) is theoretically possible.
C
C  The arguments X(*), Y(*), and Z(*), and the variable R in COMMON are
C  all INTEGER arrays of size 30.  See the comments in the routine
C  MPBLAS for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK, MPERR, MPMLP, MPNZR
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPMUL
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), Y(*), Z(*), RS, RE, XI, C, RI
C***FIRST EXECUTABLE STATEMENT  MPMUL
      CALL MPCHK (1, 4)
      I2 = T + 4
      I2P = I2 + 1
C FORM SIGN OF PRODUCT
      RS = X(1)*Y(1)
      IF (RS.NE.0) GO TO 10
C SET RESULT TO ZERO
      Z(1) = 0
      RETURN
C FORM EXPONENT OF PRODUCT
   10 RE = X(2) + Y(2)
C CLEAR ACCUMULATOR
      DO 20 I = 1, I2
   20 R(I) = 0
C PERFORM MULTIPLICATION
      C = 8
      DO 40 I = 1, T
      XI = X(I+2)
C FOR SPEED, PUT THE NUMBER WITH MANY ZEROS FIRST
      IF (XI.EQ.0) GO TO 40
      CALL MPMLP (R(I+1), Y(3), XI, MIN (T, I2 - I))
      C = C - 1
      IF (C.GT.0) GO TO 40
C CHECK FOR LEGAL BASE B DIGIT
      IF ((XI.LT.0).OR.(XI.GE.B)) GO TO 90
C PROPAGATE CARRIES AT END AND EVERY EIGHTH TIME,
C FASTER THAN DOING IT EVERY TIME.
      DO 30 J = 1, I2
      J1 = I2P - J
      RI = R(J1) + C
      IF (RI.LT.0) GO TO 70
      C = RI/B
   30 R(J1) = RI - B*C
      IF (C.NE.0) GO TO 90
      C = 8
   40 CONTINUE
      IF (C.EQ.8) GO TO 60
      IF ((XI.LT.0).OR.(XI.GE.B)) GO TO 90
      C = 0
      DO 50 J = 1, I2
      J1 = I2P - J
      RI = R(J1) + C
      IF (RI.LT.0) GO TO 70
      C = RI/B
   50 R(J1) = RI - B*C
      IF (C.NE.0) GO TO 90
C NORMALIZE AND ROUND RESULT
   60 CALL MPNZR (RS, RE, Z, 0)
      RETURN
   70 WRITE (LUN, 80)
   80 FORMAT (' *** INTEGER OVERFLOW IN MPMUL, B TOO LARGE ***')
      GO TO 110
   90 WRITE (LUN, 100)
  100 FORMAT (' *** ILLEGAL BASE B DIGIT IN CALL TO MPMUL,',
     1        ' POSSIBLE OVERWRITING PROBLEM ***')
  110 CALL MPERR
      Z(1) = 0
      RETURN
      END
