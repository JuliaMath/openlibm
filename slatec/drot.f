*DECK DROT
      SUBROUTINE DROT (N, DX, INCX, DY, INCY, DC, DS)
C***BEGIN PROLOGUE  DROT
C***PURPOSE  Apply a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A8
C***TYPE      DOUBLE PRECISION (SROT-S, DROT-D, CSROT-C)
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
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C       DC  D.P. element of rotation matrix
C       DS  D.P. element of rotation matrix
C
C     --Output--
C       DX  rotated vector DX (unchanged if N .LE. 0)
C       DY  rotated vector DY (unchanged if N .LE. 0)
C
C     Multiply the 2 x 2 matrix  ( DC DS) times the 2 x N matrix (DX**T)
C                                (-DS DC)                        (DY**T)
C     where **T indicates transpose.  The elements of DX are in
C     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY.
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
C***END PROLOGUE  DROT
      DOUBLE PRECISION DX, DY, DC, DS, ZERO, ONE, W, Z
      DIMENSION DX(*), DY(*)
      SAVE ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C***FIRST EXECUTABLE STATEMENT  DROT
      IF (N .LE. 0 .OR. (DS .EQ. ZERO .AND. DC .EQ. ONE)) GO TO 40
      IF (.NOT. (INCX .EQ. INCY .AND. INCX .GT. 0)) GO TO 20
C
C          Code for equal and positive increments.
C
           NSTEPS=INCX*N
           DO 10 I = 1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=DC*W+DS*Z
                DY(I)=-DS*W+DC*Z
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
                W=DX(KX)
                Z=DY(KY)
                DX(KX)=DC*W+DS*Z
                DY(KY)=-DS*W+DC*Z
                KX=KX+INCX
                KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
C
      RETURN
      END
