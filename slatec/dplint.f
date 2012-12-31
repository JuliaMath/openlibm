*DECK DPLINT
      SUBROUTINE DPLINT (N, X, Y, C)
C***BEGIN PROLOGUE  DPLINT
C***PURPOSE  Produce the polynomial which interpolates a set of discrete
C            data points.
C***LIBRARY   SLATEC
C***CATEGORY  E1B
C***TYPE      DOUBLE PRECISION (POLINT-S, DPLINT-D)
C***KEYWORDS  POLYNOMIAL INTERPOLATION
C***AUTHOR  Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Abstract
C        Subroutine DPLINT is designed to produce the polynomial which
C     interpolates the data  (X(I),Y(I)), I=1,...,N.  DPLINT sets up
C     information in the array C which can be used by subroutine DPOLVL
C     to evaluate the polynomial and its derivatives and by subroutine
C     DPOLCF to produce the coefficients.
C
C     Formal Parameters
C     *** All TYPE REAL variables are DOUBLE PRECISION ***
C     N  - the number of data points  (N .GE. 1)
C     X  - the array of abscissas (all of which must be distinct)
C     Y  - the array of ordinates
C     C  - an array of information used by subroutines
C     *******  Dimensioning Information  *******
C     Arrays X,Y, and C must be dimensioned at least N in the calling
C     program.
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPLINT
      INTEGER I,K,KM1,N
      DOUBLE PRECISION DIF,C(*),X(*),Y(*)
C***FIRST EXECUTABLE STATEMENT  DPLINT
      IF (N .LE. 0) GO TO 91
      C(1)=Y(1)
      IF(N .EQ. 1) RETURN
      DO 10010 K=2,N
      C(K)=Y(K)
      KM1=K-1
      DO 10010 I=1,KM1
C     CHECK FOR DISTINCT X VALUES
      DIF = X(I)-X(K)
      IF (DIF .EQ. 0.0) GO TO 92
      C(K) = (C(I)-C(K))/DIF
10010 CONTINUE
      RETURN
   91 CALL XERMSG ('SLATEC', 'DPLINT', 'N IS ZERO OR NEGATIVE.', 2, 1)
      RETURN
   92 CALL XERMSG ('SLATEC', 'DPLINT',
     +   'THE ABSCISSAS ARE NOT DISTINCT.', 2, 1)
      RETURN
      END
