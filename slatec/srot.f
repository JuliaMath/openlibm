*DECK SROT
      SUBROUTINE SROT (N, SX, INCX, SY, INCY, SC, SS)
C***BEGIN PROLOGUE  SROT
C***PURPOSE  Apply a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A8
C***TYPE      SINGLE PRECISION (SROT-S, DROT-D, CSROT-C)
C***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
C             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
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
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C       SC  element of rotation matrix
C       SS  element of rotation matrix
C
C     --Output--
C       SX  rotated vector SX (unchanged if N .LE. 0)
C       SY  rotated vector SY (unchanged if N .LE. 0)
C
C     Multiply the 2 x 2 matrix  ( SC SS) times the 2 x N matrix (SX**T)
C                                (-SS SC)                        (SY**T)
C     where **T indicates transpose.  The elements of SX are in
C     SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SROT
      REAL SX, SY, SC, SS, ZERO, ONE, W, Z
      DIMENSION SX(*), SY(*)
      SAVE ZERO, ONE
      DATA ZERO, ONE /0.0E0, 1.0E0/
C***FIRST EXECUTABLE STATEMENT  SROT
      IF (N .LE. 0 .OR. (SS .EQ. ZERO .AND. SC .EQ. ONE)) GO TO 40
      IF (.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
C          Code for equal and positive increments.
C
           NSTEPS=INCX*N
           DO 10 I = 1,NSTEPS,INCX
                W=SX(I)
                Z=SY(I)
                SX(I)=SC*W+SS*Z
                SY(I)=-SS*W+SC*Z
   10           CONTINUE
           GO TO 40
C
C     Code for unequal or nonpositive increments.
C
   20 CONTINUE
           KX=1
           KY=1
C
           IF (INCX .LT. 0) KX = 1-(N-1)*INCX
           IF (INCY .LT. 0) KY = 1-(N-1)*INCY
C
           DO 30 I = 1,N
                W=SX(KX)
                Z=SY(KY)
                SX(KX)=SC*W+SS*Z
                SY(KY)=-SS*W+SC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END
