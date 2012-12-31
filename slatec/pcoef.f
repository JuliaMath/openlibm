*DECK PCOEF
      SUBROUTINE PCOEF (L, C, TC, A)
C***BEGIN PROLOGUE  PCOEF
C***PURPOSE  Convert the POLFIT coefficients to Taylor series form.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A2
C***TYPE      SINGLE PRECISION (PCOEF-S, DPCOEF-D)
C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C***DESCRIPTION
C
C     Written BY L. F. Shampine and S. M. Davenport.
C
C     Abstract
C
C     POLFIT  computes the least squares polynomial fit of degree  L  as
C     a sum of orthogonal polynomials.  PCOEF  changes this fit to its
C     Taylor expansion about any point  C , i.e. writes the polynomial
C     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial
C     in powers of X, but a suitable non-zero  C  often leads to
C     polynomials which are better scaled and more accurately evaluated.
C
C     The parameters for  PCOEF  are
C
C     INPUT --
C         L -      Indicates the degree of polynomial to be changed to
C                  its Taylor expansion.  To obtain the Taylor
C                  coefficients in reverse order, input  L  as the
C                  negative of the degree desired.  The absolute value
C                  of L  must be less than or equal to NDEG, the highest
C                  degree polynomial fitted by  POLFIT .
C         C -      The point about which the Taylor expansion is to be
C                  made.
C         A -      Work and output array containing values from last
C                  call to  POLFIT .
C
C     OUTPUT --
C         TC -     Vector containing the first LL+1 Taylor coefficients
C                  where LL=ABS(L).  If  L.GT.0 , the coefficients are
C                  in the usual Taylor series order, i.e.
C                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N
C                  If L .LT. 0, the coefficients are in reverse order,
C                  i.e.
C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  PVALUE
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  PCOEF
C
      DIMENSION A(*), TC(*)
C***FIRST EXECUTABLE STATEMENT  PCOEF
      LL = ABS(L)
      LLP1 = LL + 1
      CALL PVALUE (LL,LL,C,TC(1),TC(2),A)
      IF (LL .LT. 2) GO TO 2
      FAC = 1.0
      DO 1 I = 3,LLP1
        FAC = FAC*(I-1)
 1      TC(I) = TC(I)/FAC
 2    IF (L .GE. 0) GO TO 4
      NR = LLP1/2
      LLP2 = LL + 2
      DO 3 I = 1,NR
        SAVE = TC(I)
        NEW = LLP2 - I
        TC(I) = TC(NEW)
 3      TC(NEW) = SAVE
 4    RETURN
      END
