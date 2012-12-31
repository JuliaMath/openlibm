*DECK SROTMG
      SUBROUTINE SROTMG (SD1, SD2, SX1, SY1, SPARAM)
C***BEGIN PROLOGUE  SROTMG
C***PURPOSE  Construct a modified Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      SINGLE PRECISION (SROTMG-S, DROTMG-D)
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
C      SD1  single precision scalar
C      SD2  single precision scalar
C      SX1  single precision scalar
C      SY2  single precision scalar
C   SPARAM  S.P. 5-vector. SPARAM(1)=SFLAG defined below.
C           Locations 2-5 contain the rotation matrix.
C
C     --Output--
C      SD1  changed to represent the effect of the transformation
C      SD2  changed to represent the effect of the transformation
C      SX1  changed to represent the effect of the transformation
C      SY2  unchanged
C
C     Construct the modified Givens transformation matrix H which zeros
C     the second component of the 2-vector  (SQRT(SD1)*SX1,SQRT(SD2)*
C     SY2)**T.
C     With SPARAM(1)=SFLAG, H has one of the following forms:
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C
C     Locations 2-5 of SPARAM contain SH11, SH21, SH12, and SH22,
C     respectively.  (Values of 1.E0, -1.E0, or 0.E0 implied by the
C     value of SPARAM(1) are not stored in SPARAM.)
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780301  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920316  Prologue corrected.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SROTMG
      DIMENSION SPARAM(5)
      SAVE ZERO, ONE, TWO, GAM, GAMSQ, RGAMSQ
      DATA ZERO, ONE, TWO /0.0E0, 1.0E0, 2.0E0/
      DATA GAM, GAMSQ, RGAMSQ /4096.0E0, 1.67772E7, 5.96046E-8/
C***FIRST EXECUTABLE STATEMENT  SROTMG
      IF (.NOT. SD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-SX1..
          GO TO 60
   10 CONTINUE
C     CASE-SD1-NONNEGATIVE
      SP2=SD2*SY1
      IF (.NOT. SP2 .EQ. ZERO) GO TO 20
          SFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      SP1=SD1*SX1
      SQ2=SP2*SY1
      SQ1=SP1*SX1
C
      IF (.NOT. ABS(SQ1) .GT. ABS(SQ2)) GO TO 40
          SH21=-SY1/SX1
          SH12=SP2/SP1
C
          SU=ONE-SH12*SH21
C
          IF (.NOT. SU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   30     CONTINUE
               SFLAG=ZERO
               SD1=SD1/SU
               SD2=SD2/SU
               SX1=SX1*SU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF (.NOT. SQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   50     CONTINUE
               SFLAG=ONE
               SH11=SP1/SP2
               SH22=SX1/SY1
               SU=ONE+SH11*SH22
               STEMP=SD2/SU
               SD2=SD1/SU
               SD1=STEMP
               SX1=SY1*SU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
          SFLAG=-ONE
          SH11=ZERO
          SH12=ZERO
          SH21=ZERO
          SH22=ZERO
C
          SD1=ZERO
          SD2=ZERO
          SX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF (.NOT. SFLAG .GE. ZERO) GO TO 90
C
          IF (.NOT. SFLAG .EQ. ZERO) GO TO 80
          SH11=ONE
          SH22=ONE
          SFLAG=-ONE
          GO TO 90
   80     CONTINUE
          SH21=-ONE
          SH12=ONE
          SFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF (.NOT. SD1 .LE. RGAMSQ) GO TO 130
               IF (SD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               SD1=SD1*GAM**2
               SX1=SX1/GAM
               SH11=SH11/GAM
               SH12=SH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF (.NOT. SD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               SD1=SD1/GAM**2
               SX1=SX1*GAM
               SH11=SH11*GAM
               SH12=SH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF (.NOT. ABS(SD2) .LE. RGAMSQ) GO TO 190
               IF (SD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               SD2=SD2*GAM**2
               SH21=SH21/GAM
               SH22=SH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF (.NOT. ABS(SD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               SD2=SD2/GAM**2
               SH21=SH21*GAM
               SH22=SH22*GAM
          GO TO 200
  220 CONTINUE
          IF (SFLAG) 250,230,240
  230     CONTINUE
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               GO TO 260
  240     CONTINUE
               SPARAM(2)=SH11
               SPARAM(5)=SH22
               GO TO 260
  250     CONTINUE
               SPARAM(2)=SH11
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               SPARAM(5)=SH22
  260 CONTINUE
          SPARAM(1)=SFLAG
          RETURN
      END
