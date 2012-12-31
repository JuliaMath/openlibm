*DECK AVINT
      SUBROUTINE AVINT (X, Y, N, XLO, XUP, ANS, IERR)
C***BEGIN PROLOGUE  AVINT
C***PURPOSE  Integrate a function tabulated at arbitrarily spaced
C            abscissas using overlapping parabolas.
C***LIBRARY   SLATEC
C***CATEGORY  H2A1B2
C***TYPE      SINGLE PRECISION (AVINT-S, DAVINT-D)
C***KEYWORDS  INTEGRATION, QUADRATURE, TABULATED DATA
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C         AVINT integrates a function tabulated at arbitrarily spaced
C         abscissas.  The limits of integration need not coincide
C         with the tabulated abscissas.
C
C         A method of overlapping parabolas fitted to the data is used
C         provided that there are at least 3 abscissas between the
C         limits of integration.  AVINT also handles two special cases.
C         If the limits of integration are equal, AVINT returns a result
C         of zero regardless of the number of tabulated values.
C         If there are only two function values, AVINT uses the
C         trapezoid rule.
C
C     Description of Parameters
C         The user must dimension all arrays appearing in the call list
C              X(N), Y(N).
C
C         Input--
C         X    - real array of abscissas, which must be in increasing
C                order.
C         Y    - real array of functional values. i.e., Y(I)=FUNC(X(I)).
C         N    - the integer number of function values supplied.
C                N .GE. 2 unless XLO = XUP.
C         XLO  - real lower limit of integration.
C         XUP  - real upper limit of integration.
C                Must have XLO .LE. XUP.
C
C         Output--
C         ANS  - computed approximate value of integral
C         IERR - a status code
C              --normal code
C                =1 means the requested integration was performed.
C              --abnormal codes
C                =2 means XUP was less than XLO.
C                =3 means the number of X(I) between XLO and XUP
C                   (inclusive) was less than 3 and neither of the two
C                   special cases described in the Abstract occurred.
C                   No integration was performed.
C                =4 means the restriction X(I+1) .GT. X(I) was violated.
C                =5 means the number N of function values was .LT. 2.
C                ANS is set to zero if IERR=2,3,4,or 5.
C
C     AVINT is documented completely in SC-M-69-335
C     Original program from "Numerical Integration" by Davis &
C     Rabinowitz.
C     Adaptation and modifications for Sandia Mathematical Program
C     Library by Rondall E. Jones.
C
C***REFERENCES  R. E. Jones, Approximate integrator of functions
C                 tabulated at arbitrarily spaced abscissas,
C                 Report SC-M-69-335, Sandia Laboratories, 1969.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   690901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  AVINT
C
      DOUBLE PRECISION R3,RP5,SUM,SYL,SYL2,SYL3,SYU,SYU2,SYU3,X1,X2,X3
     1,X12,X13,X23,TERM1,TERM2,TERM3,A,B,C,CA,CB,CC
      DIMENSION X(*),Y(*)
C***FIRST EXECUTABLE STATEMENT  AVINT
      IERR=1
      ANS =0.0
      IF (XLO-XUP) 3,100,200
    3 IF (N.LT.2) GO TO 215
      DO 5 I=2,N
      IF (X(I).LE.X(I-1)) GO TO 210
      IF (X(I).GT.XUP) GO TO 6
    5 CONTINUE
    6 CONTINUE
      IF (N.GE.3) GO TO 9
C
C     SPECIAL N=2 CASE
      SLOPE = (Y(2)-Y(1))/(X(2)-X(1))
      FL = Y(1) + SLOPE*(XLO-X(1))
      FR = Y(2) + SLOPE*(XUP-X(2))
      ANS = 0.5*(FL+FR)*(XUP-XLO)
      RETURN
    9 CONTINUE
      IF (X(N-2).LT.XLO)  GO TO 205
      IF (X(3).GT.XUP)    GO TO 205
      I = 1
   10 IF (X(I).GE.XLO) GO TO 15
      I = I+1
      GO TO 10
   15 INLFT = I
      I = N
   20 IF (X(I).LE.XUP) GO TO 25
      I = I-1
      GO TO 20
   25 INRT = I
      IF ((INRT-INLFT).LT.2) GO TO 205
      ISTART = INLFT
      IF (INLFT.EQ.1) ISTART = 2
      ISTOP  = INRT
      IF (INRT.EQ.N)  ISTOP  = N-1
C
      R3 = 3.0D0
      RP5= 0.5D0
      SUM = 0.0
      SYL = XLO
      SYL2= SYL*SYL
      SYL3= SYL2*SYL
C
      DO 50 I=ISTART,ISTOP
      X1 = X(I-1)
      X2 = X(I)
      X3 = X(I+1)
      X12 = X1-X2
      X13 = X1-X3
      X23 = X2-X3
      TERM1 = DBLE(Y(I-1))/(X12*X13)
      TERM2 =-DBLE(Y(I)) /(X12*X23)
      TERM3 = DBLE(Y(I+1))/(X13*X23)
      A = TERM1+TERM2+TERM3
      B = -(X2+X3)*TERM1 - (X1+X3)*TERM2 - (X1+X2)*TERM3
      C = X2*X3*TERM1 + X1*X3*TERM2 + X1*X2*TERM3
      IF (I-ISTART) 30,30,35
   30 CA = A
      CB = B
      CC = C
      GO TO 40
   35 CA = 0.5*(A+CA)
      CB = 0.5*(B+CB)
      CC = 0.5*(C+CC)
   40 SYU = X2
      SYU2= SYU*SYU
      SYU3= SYU2*SYU
      SUM = SUM + CA*(SYU3-SYL3)/R3  + CB*RP5*(SYU2-SYL2) + CC*(SYU-SYL)
      CA  = A
      CB  = B
      CC  = C
      SYL = SYU
      SYL2= SYU2
      SYL3= SYU3
   50 CONTINUE
      SYU = XUP
      ANS = SUM + CA*(SYU**3-SYL3)/R3 + CB*RP5*(SYU**2-SYL2)
     1  + CC*(SYU-SYL)
  100 RETURN
  200 IERR=2
      CALL XERMSG ('SLATEC', 'AVINT',
     +   'THE UPPER LIMIT OF INTEGRATION WAS NOT GREATER THAN THE ' //
     +   'LOWER LIMIT.', 4, 1)
      RETURN
  205 IERR=3
      CALL XERMSG ('SLATEC', 'AVINT',
     +   'THERE WERE LESS THAN THREE FUNCTION VALUES BETWEEN THE ' //
     +   'LIMITS OF INTEGRATION.', 4, 1)
      RETURN
  210 IERR=4
      CALL XERMSG ('SLATEC', 'AVINT',
     +   'THE ABSCISSAS WERE NOT STRICTLY INCREASING.  MUST HAVE ' //
     +   'X(I-1) .LT. X(I) FOR ALL I.', 4, 1)
      RETURN
  215 IERR=5
      CALL XERMSG ('SLATEC', 'AVINT',
     +   'LESS THAN TWO FUNCTION VALUES WERE SUPPLIED.', 4, 1)
      RETURN
      END
