*DECK DQDOTA
      DOUBLE PRECISION FUNCTION DQDOTA (N, DB, QC, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DQDOTA
C***PURPOSE  Compute the inner product of two vectors with extended
C            precision accumulation and result.
C***LIBRARY   SLATEC
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (DQDOTA-D)
C***KEYWORDS  DOT PRODUCT, INNER PRODUCT
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C               B L A S  Subprogram
C   Description of Parameters
C
C    --Input--
C       N  number of elements in input vector(S)
C      DB  double precision scalar to be added to inner product
C      QC  extended precision scalar to be added to inner product
C      DX  double precision vector with N elements
C    INCX  storage spacing between elements of DX
C      DY  double precision vector with N elements
C    INCY  storage spacing between elements of DY
C
C    --Output--
C  DQDOTA  double precision result
C      QC  extended precision result
C
C    D.P. dot product with extended precision accumulation (and result)
C    QC and DQDOTA are set = DB + QC + sum for I = 0 to N-1 of
C      DX(LX+I*INCX) * DY(LY+I*INCY),  where QC is an extended
C      precision result previously computed by DQDOTI or DQDOTA
C      and LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C      defined in a similar way using INCY.  The MP package by
C      Richard P. Brent is used for the extended precision arithmetic.
C
C    Fred T. Krogh,  JPL,  1977,  June 1
C
C    The common block for the MP package is name MPCOM.  If local
C    variable I1 is zero, DQDOTA calls MPBLAS to initialize
C    the MP package and reset I1 to 1.
C
C    The argument QC(*) and the local variables QX and QY are INTEGER
C    arrays of size 30.  See the comments in the routine MPBLAS for the
C    reason for this choice.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  MPADD, MPBLAS, MPCDM, MPCMD, MPMUL
C***COMMON BLOCKS    MPCOM
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   930124  Increased Array sizes for SUN -r8.  (RWC)
C***END PROLOGUE  DQDOTA
      DOUBLE PRECISION DX(*), DY(*), DB
      INTEGER  QC(30), QX(30), QY(30)
      COMMON /MPCOM/  MPB, MPT, MPM, MPLUN, MPMXR, MPR(30)
      SAVE I1
      DATA  I1 / 0 /
C***FIRST EXECUTABLE STATEMENT  DQDOTA
      IF (I1 .EQ. 0) CALL MPBLAS(I1)
      IF (DB .EQ. 0.D0) GO TO 20
      CALL MPCDM(DB, QX)
      CALL MPADD(QC, QX, QC)
   20 IF (N .EQ. 0) GO TO 40
      IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N + 1) * INCX + 1
      IF (INCY .LT. 0) IY = (-N + 1) * INCY + 1
      DO  30  I = 1,N
         CALL MPCDM(DX(IX), QX)
         CALL MPCDM(DY(IY), QY)
         CALL MPMUL(QX, QY, QX)
         CALL MPADD(QC, QX, QC)
         IX = IX + INCX
         IY = IY + INCY
   30 CONTINUE
   40 CALL MPCMD(QC, DQDOTA)
      RETURN
      END
