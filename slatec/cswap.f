*DECK CSWAP
      SUBROUTINE CSWAP (N, CX, INCX, CY, INCY)
C***BEGIN PROLOGUE  CSWAP
C***PURPOSE  Interchange two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      COMPLEX (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
C***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C       CY  complex vector with N elements
C     INCY  storage spacing between elements of CY
C
C     --Output--
C       CX  input vector CY (unchanged if N .LE. 0)
C       CY  input vector CX (unchanged if N .LE. 0)
C
C     Interchange complex CX and complex CY
C     For I = 0 to N-1, interchange  CX(LX+I*INCX) and CY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSWAP
      COMPLEX CX(*),CY(*),CTEMP
C***FIRST EXECUTABLE STATEMENT  CSWAP
      IF (N .LE. 0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) GO TO 20
C
C     Code for unequal or nonpositive increments.
C
      KX = 1
      KY = 1
      IF (INCX .LT. 0) KX = 1+(1-N)*INCX
      IF (INCY .LT. 0) KY = 1+(1-N)*INCY
      DO 10 I = 1,N
        CTEMP = CX(KX)
        CX(KX) = CY(KY)
        CY(KY) = CTEMP
        KX = KX + INCX
        KY = KY + INCY
   10 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   20 NS = N*INCX
      DO 30 I = 1,NS,INCX
        CTEMP = CX(I)
        CX(I) = CY(I)
        CY(I) = CTEMP
   30 CONTINUE
      RETURN
      END
