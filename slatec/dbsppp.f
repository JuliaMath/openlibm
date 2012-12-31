*DECK DBSPPP
      SUBROUTINE DBSPPP (T, A, N, K, LDC, C, XI, LXI, WORK)
C***BEGIN PROLOGUE  DBSPPP
C***PURPOSE  Convert the B-representation of a B-spline to the piecewise
C            polynomial (PP) form.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (BSPPP-S, DBSPPP-D)
C***KEYWORDS  B-SPLINE, PIECEWISE POLYNOMIAL
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract    **** a double precision routine ****
C         DBSPPP is the BSPLPP routine of the reference.
C
C         DBSPPP converts the B-representation (T,A,N,K) to the
C         piecewise polynomial (PP) form (C,XI,LXI,K) for use with
C         DPPVAL.  Here XI(*), the break point array of length LXI, is
C         the knot array T(*) with multiplicities removed.  The columns
C         of the matrix C(I,J) contain the right Taylor derivatives
C         for the polynomial expansion about XI(J) for the intervals
C         XI(J) .LE. X .LE. XI(J+1), I=1,K, J=1,LXI.  Function DPPVAL
C         makes this evaluation at a specified point X in
C         XI(1) .LE. X .LE. XI(LXI+1)
C
C     Description of Arguments
C
C         Input      T,A are double precision
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          LDC     - leading dimension of C, LDC .GE. K
C
C         Output     C,XI,WORK are double precision
C          C       - matrix of dimension at least (K,LXI) containing
C                    right derivatives at break points
C          XI      - XI break point vector of length LXI+1
C          LXI     - number of break points, LXI .LE. N-K+1
C          WORK    - work vector of length K*(N+3)
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  DBSPDR, DBSPEV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBSPPP
C
      INTEGER ILEFT, INEV, K, LDC, LXI, N, NK
      DOUBLE PRECISION A, C, T, WORK, XI
C     DIMENSION T(N+K),XI(LXI+1),C(LDC,*)
C     HERE, * = THE FINAL VALUE OF THE OUTPUT PARAMETER LXI.
      DIMENSION T(*), A(*), WORK(*), XI(*), C(LDC,*)
C***FIRST EXECUTABLE STATEMENT  DBSPPP
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      IF(LDC.LT.K) GO TO 110
      CALL DBSPDR(T, A, N, K, K, WORK)
      LXI = 0
      XI(1) = T(K)
      INEV = 1
      NK = N*K + 1
      DO 10 ILEFT=K,N
        IF (T(ILEFT+1).EQ.T(ILEFT)) GO TO 10
        LXI = LXI + 1
        XI(LXI+1) = T(ILEFT+1)
        CALL DBSPEV(T,WORK(1),N,K, K,XI(LXI),INEV,C(1,LXI),WORK(NK))
   10 CONTINUE
      RETURN
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPPP', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPPP', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBSPPP', 'LDC DOES NOT SATISFY LDC.GE.K',
     +   2, 1)
      RETURN
      END
