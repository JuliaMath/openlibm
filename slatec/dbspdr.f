*DECK DBSPDR
      SUBROUTINE DBSPDR (T, A, N, K, NDERIV, AD)
C***BEGIN PROLOGUE  DBSPDR
C***PURPOSE  Use the B-representation to construct a divided difference
C            table preparatory to a (right) derivative calculation.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (BSPDR-S, DBSPDR-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, DIFFERENTIATION OF SPLINES,
C             INTERPOLATION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract     **** a double precision routine ****
C         DBSPDR is the BSPLDR routine of the reference.
C
C         DBSPDR uses the B-representation (T,A,N,K) to construct a
C         divided difference table ADIF preparatory to a (right)
C         derivative calculation in DBSPEV.  The lower triangular matrix
C         ADIF is stored in vector AD by columns.  The arrays are
C         related by
C
C           ADIF(I,J) = AD(I-J+1 + (2*N-J+2)*(J-1)/2)
C
C         I = J,N   ,   J=1,NDERIV.
C
C     Description of Arguments
C
C         Input      T,A are double precision
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the spline, K .GE. 1
C          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K.
C                    NDERIV=1 gives the zero-th derivative =
C                    function value
C
C         Output     AD is double precision
C          AD      - table of differences in a vector of length
C                    (2*N-NDERIV+1)*NDERIV/2 for input to DBSPEV
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBSPDR
C
C
      INTEGER I, ID, II, IPKMID, JJ, JM, K, KMID, N, NDERIV
      DOUBLE PRECISION A, AD, DIFF, FKMID, T
C     DIMENSION T(N+K), AD((2*N-NDERIV+1)*NDERIV/2)
      DIMENSION T(*), A(*), AD(*)
C***FIRST EXECUTABLE STATEMENT  DBSPDR
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 110
      DO 10 I=1,N
        AD(I) = A(I)
   10 CONTINUE
      IF (NDERIV.EQ.1) RETURN
      KMID = K
      JJ = N
      JM = 0
      DO 30 ID=2,NDERIV
        KMID = KMID - 1
        FKMID = KMID
        II = 1
        DO 20 I=ID,N
          IPKMID = I + KMID
          DIFF = T(IPKMID) - T(I)
          IF (DIFF.NE.0.0D0) AD(II+JJ) = (AD(II+JM+1)-AD(II+JM))/
     1     DIFF*FKMID
          II = II + 1
   20   CONTINUE
        JM = JJ
        JJ = JJ + N - ID + 1
   30 CONTINUE
      RETURN
C
C
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPDR', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPDR', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPDR',
     +   'NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K', 2, 1)
      RETURN
      END
