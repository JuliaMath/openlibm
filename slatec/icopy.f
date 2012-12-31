*DECK ICOPY
      SUBROUTINE ICOPY (N, IX, INCX, IY, INCY)
C***BEGIN PROLOGUE  ICOPY
C***PURPOSE  Copy a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      INTEGER (ICOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Boland, W. Robert, (LANL)
C           Clemens, Reginald, (PLK)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       IX  integer vector with N elements
C     INCX  storage spacing between elements of IX
C       IY  integer vector with N elements
C     INCY  storage spacing between elements of IY
C
C     --Output--
C       IY  copy of vector IX (unchanged if N .LE. 0)
C
C     Copy integer IX to integer IY.
C     For I = 0 to N-1, copy  IX(LX+I*INCX) to IY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930201  DATE WRITTEN
C***END PROLOGUE  ICOPY
      INTEGER IX(*), IY(*)
C***FIRST EXECUTABLE STATEMENT  ICOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IIX = 1
      IIY = 1
      IF (INCX .LT. 0) IIX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IIY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        IY(IIY) = IX(IIX)
        IIX = IIX + INCX
        IIY = IIY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        IY(I) = IX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        IY(I) = IX(I)
        IY(I+1) = IX(I+1)
        IY(I+2) = IX(I+2)
        IY(I+3) = IX(I+3)
        IY(I+4) = IX(I+4)
        IY(I+5) = IX(I+5)
        IY(I+6) = IX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        IY(I) = IX(I)
   70 CONTINUE
      RETURN
      END
