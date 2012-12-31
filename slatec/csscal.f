*DECK CSSCAL
      SUBROUTINE CSSCAL (N, SA, CX, INCX)
C***BEGIN PROLOGUE  CSSCAL
C***PURPOSE  Scale a complex vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      COMPLEX (CSSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
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
C       SA  single precision scale factor
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C       CX  scaled result (unchanged if N .LE. 0)
C
C     Replace complex CX by (single precision SA) * (complex CX)
C     For I = 0 to N-1, replace CX(IX+I*INCX) with  SA * CX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
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
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSSCAL
      COMPLEX CX(*)
      REAL SA
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  CSSCAL
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CX(IX) = SA*CX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 DO 30 I = 1,N
        CX(I) = SA*CX(I)
   30 CONTINUE
      RETURN
      END
