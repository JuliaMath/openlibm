*DECK CSROT
      SUBROUTINE CSROT (N, CX, INCX, CY, INCY, C, S)
C***BEGIN PROLOGUE  CSROT
C***PURPOSE  Apply a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      COMPLEX (SROT-S, DROT-D, CSROT-C)
C***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
C             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
C***AUTHOR  Dongarra, J., (ANL)
C***DESCRIPTION
C
C     CSROT applies the complex Givens rotation
C
C          (X)   ( C S)(X)
C          (Y) = (-S C)(Y)
C
C     N times where for I = 0,...,N-1
C
C          X = CX(LX+I*INCX)
C          Y = CY(LY+I*INCY),
C
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C     Argument Description
C
C        N      (integer)  number of elements in each vector
C
C        CX     (complex array)  beginning of one vector
C
C        INCX   (integer)  memory spacing of successive elements
C               of vector CX
C
C        CY     (complex array)  beginning of the other vector
C
C        INCY   (integer)  memory spacing of successive elements
C               of vector CY
C
C        C      (real)  cosine term of the rotation
C
C        S      (real)  sine term of the rotation.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810223  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSROT
      COMPLEX CX(*), CY(*), CTEMP
      REAL C, S
      INTEGER I, INCX, INCY, IX, IY, N
C***FIRST EXECUTABLE STATEMENT  CSROT
      IF (N .LE. 0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1)GO TO 20
C
C     Code for unequal increments or equal increments not equal to 1.
C
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CTEMP = C*CX(IX) + S*CY(IY)
        CY(IY) = C*CY(IY) - S*CX(IX)
        CX(IX) = CTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
   20 DO 30 I = 1,N
        CTEMP = C*CX(I) + S*CY(I)
        CY(I) = C*CY(I) - S*CX(I)
        CX(I) = CTEMP
   30 CONTINUE
      RETURN
      END
