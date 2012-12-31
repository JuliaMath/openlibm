*DECK ISWAP
      SUBROUTINE ISWAP (N, IX, INCX, IY, INCY)
C***BEGIN PROLOGUE  ISWAP
C***PURPOSE  Interchange two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      INTEGER (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
C***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Vandevender, W. H., (SNLA)
C***DESCRIPTION
C
C                Extended B L A S  Subprogram
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
C       IX  input vector IY (unchanged if N .LE. 0)
C       IY  input vector IX (unchanged if N .LE. 0)
C
C     Interchange integer IX and integer IY.
C     For I = 0 to N-1, interchange  IX(LX+I*INCX) and IY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   850601  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  ISWAP
      INTEGER IX(*), IY(*), ITEMP1, ITEMP2, ITEMP3
C***FIRST EXECUTABLE STATEMENT  ISWAP
      IF (N .LE. 0) RETURN
      IF (INCX .NE. INCY) GO TO 5
      IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IIX = 1
      IIY = 1
      IF (INCX .LT. 0) IIX = (1-N)*INCX + 1
      IF (INCY .LT. 0) IIY = (1-N)*INCY + 1
      DO 10 I = 1,N
        ITEMP1 = IX(IIX)
        IX(IIX) = IY(IIY)
        IY(IIY) = ITEMP1
        IIX = IIX + INCX
        IIY = IIY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 3.
C
   20 M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        ITEMP1 = IX(I)
        IX(I) = IY(I)
        IY(I) = ITEMP1
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        ITEMP1 = IX(I)
        ITEMP2 = IX(I+1)
        ITEMP3 = IX(I+2)
        IX(I) = IY(I)
        IX(I+1) = IY(I+1)
        IX(I+2) = IY(I+2)
        IY(I) = ITEMP1
        IY(I+1) = ITEMP2
        IY(I+2) = ITEMP3
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        ITEMP1 = IX(I)
        IX(I) = IY(I)
        IY(I) = ITEMP1
   70 CONTINUE
      RETURN
      END
