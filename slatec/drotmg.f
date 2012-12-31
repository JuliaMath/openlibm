*DECK DROTMG
      SUBROUTINE DROTMG (DD1, DD2, DX1, DY1, DPARAM)
C***BEGIN PROLOGUE  DROTMG
C***PURPOSE  Construct a modified Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      DOUBLE PRECISION (SROTMG-S, DROTMG-D)
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
C      DD1  double precision scalar
C      DD2  double precision scalar
C      DX1  double precision scalar
C      DX2  double precision scalar
C   DPARAM  D.P. 5-vector. DPARAM(1)=DFLAG defined below.
C           Locations 2-5 contain the rotation matrix.
C
C     --Output--
C      DD1  changed to represent the effect of the transformation
C      DD2  changed to represent the effect of the transformation
C      DX1  changed to represent the effect of the transformation
C      DX2  unchanged
C
C     Construct the modified Givens transformation matrix H which zeros
C     the second component of the 2-vector  (SQRT(DD1)*DX1,SQRT(DD2)*
C     DY2)**T.
C     With DPARAM(1)=DFLAG, H has one of the following forms:
C
C     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
C
C       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
C     H=(          )    (          )    (          )    (          )
C       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
C
C     Locations 2-5 of DPARAM contain DH11, DH21, DH12, and DH22,
C     respectively.  (Values of 1.D0, -1.D0, or 0.D0 implied by the
C     value of DPARAM(1) are not stored in DPARAM.)
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920316  Prologue corrected.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DROTMG
      DOUBLE PRECISION GAM, ONE, RGAMSQ, DD1, DD2, DH11, DH12, DH21,
     1                 DH22, DPARAM, DP1, DP2, DQ1, DQ2, DU, DY1, ZERO,
     2                 GAMSQ, DFLAG, DTEMP, DX1, TWO
      DIMENSION DPARAM(5)
      SAVE ZERO, ONE, TWO, GAM, GAMSQ, RGAMSQ
      DATA ZERO, ONE, TWO /0.0D0, 1.0D0, 2.0D0/
      DATA GAM, GAMSQ, RGAMSQ /4096.0D0, 16777216.D0, 5.9604645D-8/
C***FIRST EXECUTABLE STATEMENT  DROTMG
      IF (.NOT. DD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-DX1..
          GO TO 60
   10 CONTINUE
C     CASE-DD1-NONNEGATIVE
      DP2=DD2*DY1
      IF (.NOT. DP2 .EQ. ZERO) GO TO 20
          DFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      DP1=DD1*DX1
      DQ2=DP2*DY1
      DQ1=DP1*DX1
C
      IF (.NOT. ABS(DQ1) .GT. ABS(DQ2)) GO TO 40
          DH21=-DY1/DX1
          DH12=DP2/DP1
C
          DU=ONE-DH12*DH21
C
          IF (.NOT. DU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-DX1..
               GO TO 60
   30     CONTINUE
               DFLAG=ZERO
               DD1=DD1/DU
               DD2=DD2/DU
               DX1=DX1*DU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF (.NOT. DQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-DX1..
               GO TO 60
   50     CONTINUE
               DFLAG=ONE
               DH11=DP1/DP2
               DH22=DX1/DY1
               DU=ONE+DH11*DH22
               DTEMP=DD2/DU
               DD2=DD1/DU
               DD1=DTEMP
               DX1=DY1*DU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
          DFLAG=-ONE
          DH11=ZERO
          DH12=ZERO
          DH21=ZERO
          DH22=ZERO
C
          DD1=ZERO
          DD2=ZERO
          DX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF (.NOT. DFLAG .GE. ZERO) GO TO 90
C
          IF (.NOT. DFLAG .EQ. ZERO) GO TO 80
          DH11=ONE
          DH22=ONE
          DFLAG=-ONE
          GO TO 90
   80     CONTINUE
          DH21=-ONE
          DH12=ONE
          DFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF (.NOT. DD1 .LE. RGAMSQ) GO TO 130
               IF (DD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               DD1=DD1*GAM**2
               DX1=DX1/GAM
               DH11=DH11/GAM
               DH12=DH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF (.NOT. DD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               DD1=DD1/GAM**2
               DX1=DX1*GAM
               DH11=DH11*GAM
               DH12=DH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF (.NOT. ABS(DD2) .LE. RGAMSQ) GO TO 190
               IF (DD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               DD2=DD2*GAM**2
               DH21=DH21/GAM
               DH22=DH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF (.NOT. ABS(DD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               DD2=DD2/GAM**2
               DH21=DH21*GAM
               DH22=DH22*GAM
          GO TO 200
  220 CONTINUE
          IF (DFLAG) 250,230,240
  230     CONTINUE
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               GO TO 260
  240     CONTINUE
               DPARAM(2)=DH11
               DPARAM(5)=DH22
               GO TO 260
  250     CONTINUE
               DPARAM(2)=DH11
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               DPARAM(5)=DH22
  260 CONTINUE
          DPARAM(1)=DFLAG
          RETURN
      END
