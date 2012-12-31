*DECK BSQAD
      SUBROUTINE BSQAD (T, BCOEF, N, K, X1, X2, BQUAD, WORK)
C***BEGIN PROLOGUE  BSQAD
C***PURPOSE  Compute the integral of a K-th order B-spline using the
C            B-representation.
C***LIBRARY   SLATEC
C***CATEGORY  H2A2A1, E3, K6
C***TYPE      SINGLE PRECISION (BSQAD-S, DBSQAD-D)
C***KEYWORDS  INTEGRAL OF B-SPLINES, QUADRATURE
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         BSQAD computes the integral on (X1,X2) of a K-th order
C         B-spline using the B-representation (T,BCOEF,N,K).  Orders
C         K as high as 20 are permitted by applying a 2, 6, or 10
C         point Gauss formula on subintervals of (X1,X2) which are
C         formed by included (distinct) knots.
C
C         If orders K greater than 20 are needed, use BFQAD with
C         F(X) = 1.
C
C     Description of Arguments
C         Input
C           T      - knot array of length N+K
C           BCOEF  - B-spline coefficient array of length N
C           N      - length of coefficient array
C           K      - order of B-spline, 1 .LE. K .LE. 20
C           X1,X2  - end points of quadrature interval in
C                    T(K) .LE. X .LE. T(N+1)
C
C         Output
C           BQUAD  - integral of the B-spline over (X1,X2)
C           WORK   - work vector of length 3*K
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
C                 B-splines, Report SAND79-1825, Sandia Laboratories,
C                 December 1979.
C***ROUTINES CALLED  BVALU, INTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BSQAD
C
      INTEGER I,IL1,IL2,ILO,INBV, JF,K,LEFT,M,MF,MFLAG,N, NPK, NP1
      REAL A, AA, B, BB, BCOEF, BMA, BPA, BQUAD, C1, GPTS, GWTS, GX, Q,
     1 SUM, T, TA, TB, WORK, X1, X2, Y1, Y2
      REAL BVALU
      DIMENSION T(*), BCOEF(*), GPTS(9), GWTS(9), SUM(5), WORK(*)
C
      SAVE GPTS, GWTS
      DATA GPTS(1), GPTS(2), GPTS(3), GPTS(4), GPTS(5), GPTS(6),
     1     GPTS(7), GPTS(8), GPTS(9)/
     2     5.77350269189625764E-01,     2.38619186083196909E-01,
     3     6.61209386466264514E-01,     9.32469514203152028E-01,
     4     1.48874338981631211E-01,     4.33395394129247191E-01,
     5     6.79409568299024406E-01,     8.65063366688984511E-01,
     6     9.73906528517171720E-01/
      DATA GWTS(1), GWTS(2), GWTS(3), GWTS(4), GWTS(5), GWTS(6),
     1     GWTS(7), GWTS(8), GWTS(9)/
     2     1.00000000000000000E+00,     4.67913934572691047E-01,
     3     3.60761573048138608E-01,     1.71324492379170345E-01,
     4     2.95524224714752870E-01,     2.69266719309996355E-01,
     5     2.19086362515982044E-01,     1.49451349150580593E-01,
     6     6.66713443086881376E-02/
C
C***FIRST EXECUTABLE STATEMENT  BSQAD
      BQUAD = 0.0E0
      IF(K.LT.1 .OR. K.GT.20) GO TO 65
      IF(N.LT.K) GO TO 70
      AA = MIN(X1,X2)
      BB = MAX(X1,X2)
      IF (AA.LT.T(K)) GO TO 60
      NP1 = N + 1
      IF (BB.GT.T(NP1)) GO TO 60
      IF (AA.EQ.BB) RETURN
      NPK = N + K
C     SELECTION OF 2, 6, OR 10 POINT GAUSS FORMULA
      JF = 0
      MF = 1
      IF (K.LE.4) GO TO 10
      JF = 1
      MF = 3
      IF (K.LE.12) GO TO 10
      JF = 4
      MF = 5
   10 CONTINUE
C
      DO 20 I=1,MF
        SUM(I) = 0.0E0
   20 CONTINUE
      ILO = 1
      INBV = 1
      CALL INTRV(T, NPK, AA, ILO, IL1, MFLAG)
      CALL INTRV(T, NPK, BB, ILO, IL2, MFLAG)
      IF (IL2.GE.NP1) IL2 = N
      DO 40 LEFT=IL1,IL2
        TA = T(LEFT)
        TB = T(LEFT+1)
        IF (TA.EQ.TB) GO TO 40
        A = MAX(AA,TA)
        B = MIN(BB,TB)
        BMA = 0.5E0*(B-A)
        BPA = 0.5E0*(B+A)
        DO 30 M=1,MF
          C1 = BMA*GPTS(JF+M)
          GX = -C1 + BPA
          Y2 = BVALU(T,BCOEF,N,K,0,GX,INBV,WORK)
          GX = C1 + BPA
          Y1 = BVALU(T,BCOEF,N,K,0,GX,INBV,WORK)
          SUM(M) = SUM(M) + (Y1+Y2)*BMA
   30   CONTINUE
   40 CONTINUE
      Q = 0.0E0
      DO 50 M=1,MF
        Q = Q + GWTS(JF+M)*SUM(M)
   50 CONTINUE
      IF (X1.GT.X2) Q = -Q
      BQUAD = Q
      RETURN
C
C
   60 CONTINUE
      CALL XERMSG ('SLATEC', 'BSQAD',
     +   'X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE.T(N+1)', 2, 1)
      RETURN
   65 CONTINUE
      CALL XERMSG ('SLATEC', 'BSQAD', 'K DOES NOT SATISFY 1.LE.K.LE.20'
     +   , 2, 1)
      RETURN
   70 CONTINUE
      CALL XERMSG ('SLATEC', 'BSQAD', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
      END
