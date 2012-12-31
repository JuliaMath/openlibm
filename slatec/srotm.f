*DECK SROTM
      SUBROUTINE SROTM (N, SX, INCX, SY, INCY, SPARAM)
C***BEGIN PROLOGUE  SROTM
C***PURPOSE  Apply a modified Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A8
C***TYPE      SINGLE PRECISION (SROTM-S, DROTM-D)
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
C       SX  single precision vector with N elements
C     INCX  storage spacing between elements of SX
C       SY  single precision vector with N elements
C     INCY  storage spacing between elements of SY
C   SPARAM  5-element vector. SPARAM(1) is SFLAG described below.
C           Locations 2-5 of SPARAM contain elements of the
C           transformation matrix H described below.
C
C     --Output--
C       SX  rotated vector (unchanged if N .LE. 0)
C       SY  rotated vector (unchanged if N .LE. 0)
C
C     Apply the modified Givens transformation, H, to the 2 by N matrix
C     (SX**T)
C     (SY**T) , where **T indicates transpose.  The elements of SX are
C     in SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else
C     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
C
C     With SPARAM(1)=SFLAG, H has one of the following forms:
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C
C     See SROTMG for a description of data storage in SPARAM.
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
C***END PROLOGUE  SROTM
      DIMENSION SX(*), SY(*), SPARAM(5)
      SAVE ZERO, TWO
      DATA ZERO, TWO /0.0E0, 2.0E0/
C***FIRST EXECUTABLE STATEMENT  SROTM
      SFLAG=SPARAM(1)
      IF (N.LE.0 .OR. (SFLAG+TWO.EQ.ZERO)) GO TO 140
          IF (.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
C
               NSTEPS=N*INCX
               IF (SFLAG) 50,10,30
   10          CONTINUE
               SH12=SPARAM(4)
               SH21=SPARAM(3)
                    DO 20 I = 1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W+Z*SH12
                    SY(I)=W*SH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               SH11=SPARAM(2)
               SH22=SPARAM(5)
                    DO 40 I = 1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z
                    SY(I)=-W+SH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               SH11=SPARAM(2)
               SH12=SPARAM(4)
               SH21=SPARAM(3)
               SH22=SPARAM(5)
                    DO 60 I = 1,NSTEPS,INCX
                    W=SX(I)
                    Z=SY(I)
                    SX(I)=W*SH11+Z*SH12
                    SY(I)=W*SH21+Z*SH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF (INCX .LT. 0) KX = 1+(1-N)*INCX
          IF (INCY .LT. 0) KY = 1+(1-N)*INCY
C
          IF (SFLAG) 120,80,100
   80     CONTINUE
          SH12=SPARAM(4)
          SH21=SPARAM(3)
               DO 90 I = 1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W+Z*SH12
               SY(KY)=W*SH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          SH11=SPARAM(2)
          SH22=SPARAM(5)
               DO 110 I = 1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z
               SY(KY)=-W+SH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          SH11=SPARAM(2)
          SH12=SPARAM(4)
          SH21=SPARAM(3)
          SH22=SPARAM(5)
               DO 130 I = 1,N
               W=SX(KX)
               Z=SY(KY)
               SX(KX)=W*SH11+Z*SH12
               SY(KY)=W*SH21+Z*SH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
      END
