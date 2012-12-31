*DECK MPCDM
      SUBROUTINE MPCDM (DX, Z)
C***BEGIN PROLOGUE  MPCDM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPCDM-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C Converts double-precision number DX to multiple-precision Z.
C Some numbers will not convert exactly on machines with base
C other than two, four or sixteen. This routine is not called
C by any other routine in 'mp', so may be omitted if double-
C precision is not available.
C
C The argument Z(*) and the variable R in COMMON are both INTEGER
C arrays of size 30.  See the comments in the routine MPBLAS for the
C for the reason for this choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***ROUTINES CALLED  MPCHK, MPDIVI, MPMULI, MPNZR
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPCDM
      DOUBLE PRECISION DB, DJ, DX
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, Z(*), RS, RE, TP
C***FIRST EXECUTABLE STATEMENT  MPCDM
      CALL MPCHK (1, 4)
      I2 = T + 4
C CHECK SIGN
      IF (DX) 20, 10, 30
C IF DX = 0D0 RETURN 0
   10 Z(1) = 0
      RETURN
C DX .LT. 0D0
   20 RS = -1
      DJ = -DX
      GO TO 40
C DX .GT. 0D0
   30 RS = 1
      DJ = DX
   40 IE = 0
   50 IF (DJ.LT.1D0) GO TO 60
C INCREASE IE AND DIVIDE DJ BY 16.
      IE = IE + 1
      DJ = 0.0625D0*DJ
      GO TO 50
   60 IF (DJ.GE.0.0625D0) GO TO 70
      IE = IE - 1
      DJ = 16D0*DJ
      GO TO 60
C NOW DJ IS DY DIVIDED BY SUITABLE POWER OF 16
C SET EXPONENT TO 0
   70 RE = 0
      DB = DBLE(B)
C CONVERSION LOOP (ASSUME DOUBLE-PRECISION OPS. EXACT)
      DO 80 I = 1, I2
      DJ = DB*DJ
      R(I) = INT(DJ)
   80 DJ = DJ - DBLE(R(I))
C NORMALIZE RESULT
      CALL MPNZR (RS, RE, Z, 0)
      IB = MAX(7*B*B, 32767)/16
      TP = 1
C NOW MULTIPLY BY 16**IE
      IF (IE) 90, 130, 110
   90 K = -IE
      DO 100 I = 1, K
      TP = 16*TP
      IF ((TP.LE.IB).AND.(TP.NE.B).AND.(I.LT.K)) GO TO 100
      CALL MPDIVI (Z, TP, Z)
      TP = 1
  100 CONTINUE
      RETURN
  110 DO 120 I = 1, IE
      TP = 16*TP
      IF ((TP.LE.IB).AND.(TP.NE.B).AND.(I.LT.IE)) GO TO 120
      CALL MPMULI (Z, TP, Z)
      TP = 1
  120 CONTINUE
  130 RETURN
      END
