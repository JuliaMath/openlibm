*DECK DPPQAD
      SUBROUTINE DPPQAD (LDC, C, XI, LXI, K, X1, X2, PQUAD)
C***BEGIN PROLOGUE  DPPQAD
C***PURPOSE  Compute the integral on (X1,X2) of a K-th order B-spline
C            using the piecewise polynomial (PP) representation.
C***LIBRARY   SLATEC
C***CATEGORY  H2A2A1, E3, K6
C***TYPE      DOUBLE PRECISION (PPQAD-S, DPPQAD-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract    **** a double precision routine ****
C         DPPQAD computes the integral on (X1,X2) of a K-th order
C         B-spline using the piecewise polynomial representation
C         (C,XI,LXI,K).  Here the Taylor expansion about the left
C         end point XI(J) of the J-th interval is integrated and
C         evaluated on subintervals of (X1,X2) which are formed by
C         included break points.  Integration outside (XI(1),XI(LXI+1))
C         is permitted.
C
C     Description of Arguments
C         Input      C,XI,X1,X2 are double precision
C           LDC    - leading dimension of matrix C, LDC .GE. K
C           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI
C           XI(*)  - break point array of length LXI+1
C           LXI    - number of polynomial pieces
C           K      - order of B-spline, K .GE. 1
C           X1,X2  - end points of quadrature interval, normally in
C                    XI(1) .LE. X .LE. XI(LXI+1)
C
C         Output     PQUAD is double precision
C           PQUAD  - integral of the PP representation over (X1,X2)
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
C                 B-splines, Report SAND79-1825, Sandia Laboratories,
C                 December 1979.
C***ROUTINES CALLED  DINTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPPQAD
C
      INTEGER I, II, IL, ILO, IL1, IL2, IM, K, LDC, LEFT, LXI, MF1, MF2
      DOUBLE PRECISION A,AA,BB,C,DX,FLK,PQUAD,Q,S,SS,TA,TB,X,XI,X1,X2
      DIMENSION XI(*), C(LDC,*), SS(2)
C
C***FIRST EXECUTABLE STATEMENT  DPPQAD
      PQUAD = 0.0D0
      IF(K.LT.1) GO TO 100
      IF(LXI.LT.1) GO TO 105
      IF(LDC.LT.K) GO TO 110
      AA = MIN(X1,X2)
      BB = MAX(X1,X2)
      IF (AA.EQ.BB) RETURN
      ILO = 1
      CALL DINTRV(XI, LXI, AA, ILO, IL1, MF1)
      CALL DINTRV(XI, LXI, BB, ILO, IL2, MF2)
      Q = 0.0D0
      DO 40 LEFT=IL1,IL2
        TA = XI(LEFT)
        A = MAX(AA,TA)
        IF (LEFT.EQ.1) A = AA
        TB = BB
        IF (LEFT.LT.LXI) TB = XI(LEFT+1)
        X = MIN(BB,TB)
        DO 30 II=1,2
          SS(II) = 0.0D0
          DX = X - XI(LEFT)
          IF (DX.EQ.0.0D0) GO TO 20
          S = C(K,LEFT)
          FLK = K
          IM = K - 1
          IL = IM
          DO 10 I=1,IL
            S = S*DX/FLK + C(IM,LEFT)
            IM = IM - 1
            FLK = FLK - 1.0D0
   10     CONTINUE
          SS(II) = S*DX
   20     CONTINUE
          X = A
   30   CONTINUE
        Q = Q + (SS(1)-SS(2))
   40 CONTINUE
      IF (X1.GT.X2) Q = -Q
      PQUAD = Q
      RETURN
C
C
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPQAD', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPQAD', 'LXI DOES NOT SATISFY LXI.GE.1',
     +   2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPQAD', 'LDC DOES NOT SATISFY LDC.GE.K',
     +   2, 1)
      RETURN
      END
