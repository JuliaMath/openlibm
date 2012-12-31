*DECK SSCAL
      SUBROUTINE SSCAL (N, SA, SX, INCX)
C***BEGIN PROLOGUE  SSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      SINGLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
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
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C
C     --Output--
C       SX  single precision result (unchanged if N .LE. 0)
C
C     Replace single precision SX by single precision SA*SX.
C     For I = 0 to N-1, replace SX(IX+I*INCX) with  SA * SX(IX+I*INCX),
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
C***END PROLOGUE  SSCAL
      REAL SA, SX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  SSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        SX(IX) = SA*SX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I+1) = SA*SX(I+1)
        SX(I+2) = SA*SX(I+2)
        SX(I+3) = SA*SX(I+3)
        SX(I+4) = SA*SX(I+4)
   50 CONTINUE
      RETURN
      END
