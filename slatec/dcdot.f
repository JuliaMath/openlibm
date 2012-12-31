*DECK DCDOT
      SUBROUTINE DCDOT (N, FM, CX, INCX, CY, INCY, DCR, DCI)
C***BEGIN PROLOGUE  DCDOT
C***PURPOSE  Compute the inner product of two vectors with extended
C            precision accumulation and result.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      COMPLEX (DSDOT-D, DCDOT-C)
C***KEYWORDS  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT,
C             LINEAR ALGEBRA, VECTOR
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C    Compute the dot product of 2 complex vectors, CX and CY, e.g.
C    CX DOT CY, or, CXconjugate DOT CY.  The real and imaginary
C    parts of CX and CY are converted to double precision, the dot
C    product accumulation is done in double precision and the output
C    is given as 2 double precision numbers, corresponding to the real
C    and imaginary part of the result.
C     Input
C      N:  Number of complex components of CX and CY.
C      FM: =+1.0   compute CX DOT CY.
C          =-1.0   compute CXconjugate DOT CY.
C      CX(N):
C      CY(N):  Complex arrays of length N.
C      INCX:(Integer)   Spacing of elements of CX to use
C      INCY:(Integer)   Spacing of elements of CY to use.
C     Output
C      DCR:(Double Precision) Real part of dot product.
C      DCI:(Double Precision) Imaginary part of dot product.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DCDOT
      INTEGER I, INCX, INCY, KX, KY, N
      COMPLEX CX(*), CY(*)
      DOUBLE PRECISION DCR, DCI, DT1, DT2, DT3, DT4, FM
C***FIRST EXECUTABLE STATEMENT  DCDOT
      DCR = 0.0D0
      DCI = 0.0D0
      IF (N .LE. 0) GO TO 20
C
      KX = 1
      KY = 1
      IF (INCX .LT. 0) KX = 1+(1-N)*INCX
      IF (INCY .LT. 0) KY = 1+(1-N)*INCY
      DO 10 I = 1,N
        DT1 = DBLE(REAL(CX(KX)))
        DT2 = DBLE(REAL(CY(KY)))
        DT3 = DBLE(AIMAG(CX(KX)))
        DT4 = DBLE(AIMAG(CY(KY)))
        DCR = DCR+(DT1*DT2)-FM*(DT3*DT4)
        DCI = DCI+(DT1*DT4)+FM*(DT3*DT2)
        KX = KX+INCX
        KY = KY+INCY
   10 CONTINUE
   20 RETURN
      END
