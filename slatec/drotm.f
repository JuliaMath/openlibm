*DECK DROTM
      SUBROUTINE DROTM (N, DX, INCX, DY, INCY, DPARAM)
C***BEGIN PROLOGUE  DROTM
C***PURPOSE  Apply a modified Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A8
C***TYPE      DOUBLE PRECISION (SROTM-S, DROTM-D)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
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
C   DPARAM  5-element D.P. vector.  DPARAM(1) is DFLAG described below.
C           Locations 2-5 of SPARAM contain elements of the
C           transformation matrix H described below.
C
C     --Output--
C       DX  rotated vector (unchanged if N .LE. 0)
C       DY  rotated vector (unchanged if N .LE. 0)
C
C     Apply the modified Givens transformation, H, to the 2 by N matrix
C     (DX**T)
C     (DY**T) , where **T indicates transpose.  The elements of DX are
C     in DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY.
C
C     With DPARAM(1)=DFLAG, H has one of the following forms:
C
C     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
C
C       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
C     H=(          )    (          )    (          )    (          )
C       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
C
C     See DROTMG for a description of data storage in DPARAM.
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
C***END PROLOGUE  DROTM
      DOUBLE PRECISION DFLAG, DH12, DH22, DX, TWO, Z, DH11, DH21,
     1                 DPARAM, DY, W, ZERO
      DIMENSION DX(*), DY(*), DPARAM(5)
      SAVE ZERO, TWO
      DATA ZERO, TWO /0.0D0, 2.0D0/
C***FIRST EXECUTABLE STATEMENT  DROTM
      DFLAG=DPARAM(1)
      IF (N.LE.0 .OR. (DFLAG+TWO.EQ.ZERO)) GO TO 140
          IF (.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF (DFLAG) 50,10,30
   10          CONTINUE
               DH12=DPARAM(4)
               DH21=DPARAM(3)
                    DO 20 I = 1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W+Z*DH12
                    DY(I)=W*DH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               DH11=DPARAM(2)
               DH22=DPARAM(5)
                    DO 40 I = 1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z
                    DY(I)=-W+DH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               DH11=DPARAM(2)
               DH12=DPARAM(4)
               DH21=DPARAM(3)
               DH22=DPARAM(5)
                    DO 60 I = 1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z*DH12
                    DY(I)=W*DH21+Z*DH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF (INCX .LT. 0) KX = 1+(1-N)*INCX
          IF (INCY .LT. 0) KY = 1+(1-N)*INCY
C
          IF (DFLAG) 120,80,100
   80     CONTINUE
          DH12=DPARAM(4)
          DH21=DPARAM(3)
               DO 90 I = 1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W+Z*DH12
               DY(KY)=W*DH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          DH11=DPARAM(2)
          DH22=DPARAM(5)
               DO 110 I = 1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z
               DY(KY)=-W+DH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          DH11=DPARAM(2)
          DH12=DPARAM(4)
          DH21=DPARAM(3)
          DH22=DPARAM(5)
               DO 130 I = 1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z*DH12
               DY(KY)=W*DH21+Z*DH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
      END
