*DECK CDCDOT
      COMPLEX FUNCTION CDCDOT (N, CB, CX, INCX, CY, INCY)
C***BEGIN PROLOGUE  CDCDOT
C***PURPOSE  Compute the inner product of two vectors with extended
C            precision accumulation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      COMPLEX (SDSDOT-S, CDCDOT-C)
C***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
C       CB  complex scalar to be added to inner product
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C       CY  complex vector with N elements
C     INCY  storage spacing between elements of CY
C
C     --Output--
C   CDCDOT  complex dot product (CB if N .LE. 0)
C
C     Returns complex result with dot product accumulated in D.P.
C     CDCDOT = CB + sum for I = 0 to N-1 of CX(LX+I*INCY)*CY(LY+I*INCY)
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
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CDCDOT
      INTEGER N, INCX, INCY, I, KX, KY
      COMPLEX CX(*), CY(*), CB
      DOUBLE PRECISION DSDOTR, DSDOTI, DT1, DT2, DT3, DT4
C***FIRST EXECUTABLE STATEMENT  CDCDOT
      DSDOTR = DBLE(REAL(CB))
      DSDOTI = DBLE(AIMAG(CB))
      IF (N .LE. 0) GO TO 10
      KX = 1
      KY = 1
      IF(INCX.LT.0) KX = 1+(1-N)*INCX
      IF(INCY.LT.0) KY = 1+(1-N)*INCY
      DO 5 I = 1,N
        DT1 = DBLE(REAL(CX(KX)))
        DT2 = DBLE(REAL(CY(KY)))
        DT3 = DBLE(AIMAG(CX(KX)))
        DT4 = DBLE(AIMAG(CY(KY)))
        DSDOTR = DSDOTR+(DT1*DT2)-(DT3*DT4)
        DSDOTI = DSDOTI+(DT1*DT4)+(DT3*DT2)
        KX = KX+INCX
        KY = KY+INCY
    5 CONTINUE
   10 CDCDOT = CMPLX(REAL(DSDOTR),REAL(DSDOTI))
      RETURN
      END
