*DECK DPFQAD
      SUBROUTINE DPFQAD (F, LDC, C, XI, LXI, K, ID, X1, X2, TOL, QUAD,
     +   IERR)
C***BEGIN PROLOGUE  DPFQAD
C***PURPOSE  Compute the integral on (X1,X2) of a product of a
C            function F and the ID-th derivative of a B-spline,
C            (PP-representation).
C***LIBRARY   SLATEC
C***CATEGORY  H2A2A1, E3, K6
C***TYPE      DOUBLE PRECISION (PFQAD-S, DPFQAD-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract    **** a double precision routine ****
C         DPFQAD computes the integral on (X1,X2) of a product of a
C         function F and the ID-th derivative of a B-spline, using the
C         PP-representation (C,XI,LXI,K).  (X1,X2) is normally a sub-
C         interval of XI(1) .LE. X .LE. XI(LXI+1).  An integration
C         routine, DPPGQ8 (a modification of GAUS8), integrates the
C         product on subintervals of (X1,X2) formed by the included
C         break points.  Integration outside of (XI(1),XI(LXI+1)) is
C         permitted provided F is defined.
C
C         The maximum number of significant digits obtainable in
C         DBSQAD is the smaller of 18 and the number of digits
C         carried in double precision arithmetic.
C
C     Description of arguments
C         Input      F,C,XI,X1,X2,TOL are double precision
C           F      - external function of one argument for the
C                    integrand PF(X)=F(X)*DPPVAL(LDC,C,XI,LXI,K,ID,X,
C                    INPPV)
C           LDC    - leading dimension of matrix C, LDC .GE. K
C           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI
C           XI(*)  - break point array of length LXI+1
C           LXI    - number of polynomial pieces
C           K      - order of B-spline, K .GE. 1
C           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1
C                    ID=0 gives the spline function
C           X1,X2  - end points of quadrature interval, normally in
C                    XI(1) .LE. X .LE. XI(LXI+1)
C           TOL    - desired accuracy for the quadrature, suggest
C                    10.*DTOL .LT. TOL .LE. 0.1 where DTOL is the
C                    maximum of 1.0D-18 and double precision unit
C                    roundoff for the machine = D1MACH(4)
C
C         Output     QUAD is double precision
C           QUAD   - integral of PF(X) on (X1,X2)
C           IERR   - a status code
C                    IERR=1 normal return
C                         2 some quadrature does not meet the
C                           requested tolerance
C
C     Error Conditions
C         Improper input is a fatal error.
C         Some quadrature does not meet the requested tolerance.
C
C***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
C                 B-splines, Report SAND79-1825, Sandia Laboratories,
C                 December 1979.
C***ROUTINES CALLED  D1MACH, DINTRV, DPPGQ8, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPFQAD
C
      INTEGER ID,IERR,IFLG,ILO,IL1,IL2,INPPV,K,LDC,LEFT,LXI,MF1,MF2
      DOUBLE PRECISION A,AA,ANS,B,BB,C,Q,QUAD,TA,TB,TOL,WTOL,XI,X1,X2
      DOUBLE PRECISION D1MACH, F
      DIMENSION XI(*), C(LDC,*)
      EXTERNAL F
C
C***FIRST EXECUTABLE STATEMENT  DPFQAD
      IERR = 1
      QUAD = 0.0D0
      IF(K.LT.1) GO TO 100
      IF(LDC.LT.K) GO TO 105
      IF(ID.LT.0 .OR. ID.GE.K) GO TO 110
      IF(LXI.LT.1) GO TO 115
      WTOL = D1MACH(4)
      WTOL = MAX(WTOL,1.0D-18)
      IF (TOL.LT.WTOL .OR. TOL.GT.0.1D0) GO TO 20
      AA = MIN(X1,X2)
      BB = MAX(X1,X2)
      IF (AA.EQ.BB) RETURN
      ILO = 1
      CALL DINTRV(XI, LXI, AA, ILO, IL1, MF1)
      CALL DINTRV(XI, LXI, BB, ILO, IL2, MF2)
      Q = 0.0D0
      INPPV = 1
      DO 10 LEFT=IL1,IL2
        TA = XI(LEFT)
        A = MAX(AA,TA)
        IF (LEFT.EQ.1) A = AA
        TB = BB
        IF (LEFT.LT.LXI) TB = XI(LEFT+1)
        B = MIN(BB,TB)
        CALL DPPGQ8(F,LDC,C,XI,LXI,K,ID,A,B,INPPV,TOL,ANS,IFLG)
        IF (IFLG.GT.1) IERR = 2
        Q = Q + ANS
   10 CONTINUE
      IF (X1.GT.X2) Q = -Q
      QUAD = Q
      RETURN
C
   20 CONTINUE
      CALL XERMSG ('SLATEC', 'DPFQAD',
     +   'TOL IS LESS DTOL OR GREATER THAN 0.1', 2, 1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'DPFQAD', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'DPFQAD', 'LDC DOES NOT SATISFY LDC.GE.K',
     +   2, 1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DPFQAD',
     +   'ID DOES NOT SATISFY 0.LE.ID.LT.K', 2, 1)
      RETURN
  115 CONTINUE
      CALL XERMSG ('SLATEC', 'DPFQAD', 'LXI DOES NOT SATISFY LXI.GE.1',
     +   2, 1)
      RETURN
      END
