*DECK SAXPY
      SUBROUTINE SAXPY (N, SA, SX, INCX, SY, INCY)
C***BEGIN PROLOGUE  SAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      SINGLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
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
C       SA  single precision scalar multiplier
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C
C     --Output--
C       SY  single precision result (unchanged if N .LE. 0)
C
C     Overwrite single precision SY with single precision SA*SX +SY.
C     For I = 0 to N-1, replace  SY(LY+I*INCY) with SA*SX(LX+I*INCX) +
C       SY(LY+I*INCY),
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
C***END PROLOGUE  SAXPY
      REAL SX(*), SY(*), SA
C***FIRST EXECUTABLE STATEMENT  SAXPY
      IF (N.LE.0 .OR. SA.EQ.0.0E0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I+1) = SY(I+1) + SA*SX(I+1)
        SY(I+2) = SY(I+2) + SA*SX(I+2)
        SY(I+3) = SY(I+3) + SA*SX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        SY(I) = SA*SX(I) + SY(I)
   70 CONTINUE
      RETURN
      END
