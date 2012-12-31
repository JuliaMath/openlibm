*DECK MPADD2
      SUBROUTINE MPADD2 (X, Y, Z, Y1, TRUNC)
C***BEGIN PROLOGUE  MPADD2
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DQDOTA and DQDOTI
C***LIBRARY   SLATEC
C***TYPE      ALL (MPADD2-A)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C  Called by MPADD, MPSUB etc.
C  X, Y and Z are MP numbers, Y1 and TRUNC are integers.
C  To force call by reference rather than value/result, Y1 is
C  declared as an array, but only Y1(1) is ever used.
C  Sets Z = X + Y1(1)*ABS(Y), where Y1(1) = +- Y(1).
C  If TRUNC .EQ. 0, R*-rounding is used;  otherwise, truncation.
C  R*-rounding is defined in the Kuki and Cody reference.
C
C  The arguments X(*), Y(*), and Z(*) are all INTEGER arrays of size
C  30.  See the comments in the routine MPBLAS for the reason for this
C  choice.
C
C***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
C***REFERENCES  H. Kuki and W. J. Cody, A statistical study of floating
C                 point number systems, Communications of the ACM 16, 4
C                 (April 1973), pp. 223-230.
C               R. P. Brent, On the precision attainable with various
C                 floating-point number systems, IEEE Transactions on
C                 Computers C-22, 6 (June 1973), pp. 601-607.
C               R. P. Brent, A Fortran multiple-precision arithmetic
C                 package, ACM Transactions on Mathematical Software 4,
C                 1 (March 1978), pp. 57-70.
C               R. P. Brent, MP, a Fortran multiple-precision arithmetic
C                 package, Algorithm 524, ACM Transactions on Mathema-
C                 tical Software 4, 1 (March 1978), pp. 71-81.
C***ROUTINES CALLED  MPADD3, MPCHK, MPERR, MPNZR, MPSTR
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   ??????  Modified for use with BLAS.  Blank COMMON changed to named
C           COMMON.  R given dimension 12.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C   920528  Added a REFERENCES section revised.  (WRB)
C   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
C***END PROLOGUE  MPADD2
      COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
      INTEGER B, T, R, X(*), Y(*), Z(*), Y1(*), TRUNC
      INTEGER S, ED, RS, RE
C***FIRST EXECUTABLE STATEMENT  MPADD2
      IF (X(1).NE.0) GO TO 20
   10 CALL MPSTR(Y, Z)
      Z(1) = Y1(1)
      RETURN
   20 IF (Y1(1).NE.0) GO TO 40
   30 CALL MPSTR (X, Z)
      RETURN
C COMPARE SIGNS
   40 S = X(1)*Y1(1)
      IF (ABS(S).LE.1) GO TO 60
      CALL MPCHK (1, 4)
      WRITE (LUN, 50)
   50 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN CALL TO MPADD2,',
     1        ' POSSIBLE OVERWRITING PROBLEM ***')
      CALL MPERR
      Z(1) = 0
      RETURN
C COMPARE EXPONENTS
   60 ED = X(2) - Y(2)
      MED = ABS(ED)
      IF (ED) 90, 70, 120
C EXPONENTS EQUAL SO COMPARE SIGNS, THEN FRACTIONS IF NEC.
   70 IF (S.GT.0) GO TO 100
      DO 80 J = 1, T
      IF (X(J+2) - Y(J+2)) 100, 80, 130
   80 CONTINUE
C RESULT IS ZERO
      Z(1) = 0
      RETURN
C HERE EXPONENT(Y) .GE. EXPONENT(X)
   90 IF (MED.GT.T) GO TO 10
  100 RS = Y1(1)
      RE = Y(2)
      CALL MPADD3 (X, Y, S, MED, RE)
C NORMALIZE, ROUND OR TRUNCATE, AND RETURN
  110 CALL MPNZR (RS, RE, Z, TRUNC)
      RETURN
C ABS(X) .GT. ABS(Y)
  120 IF (MED.GT.T) GO TO 30
  130 RS = X(1)
      RE = X(2)
      CALL MPADD3 (Y, X, S, MED, RE)
      GO TO 110
      END
