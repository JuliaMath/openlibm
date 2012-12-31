*DECK SROTG
      SUBROUTINE SROTG (SA, SB, SC, SS)
C***BEGIN PROLOGUE  SROTG
C***PURPOSE  Construct a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      SINGLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
C***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
C             LINEAR ALGEBRA, VECTOR
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
C       SA  single precision scalar
C       SB  single precision scalar
C
C     --Output--
C       SA  single precision result R
C       SB  single precision result Z
C       SC  single precision result
C       SS  single precision result
C
C     Construct the Givens transformation
C
C         ( SC  SS )
C     G = (        ) ,    SC**2 + SS**2 = 1 ,
C         (-SS  SC )
C
C     which zeros the second entry of the 2-vector  (SA,SB)**T.
C
C     The quantity R = (+/-)SQRT(SA**2 + SB**2) overwrites SA in
C     storage.  The value of SB is overwritten by a value Z which
C     allows SC and SS to be recovered by the following algorithm:
C
C           If Z=1  set  SC=0.0  and  SS=1.0
C           If ABS(Z) .LT. 1  set  SC=SQRT(1-Z**2)  and  SS=Z
C           If ABS(Z) .GT. 1  set  SC=1/Z  and  SS=SQRT(1-SC**2)
C
C     Normally, the subprogram SROT(N,SX,INCX,SY,INCY,SC,SS) will
C     next be called to apply the transformation to a 2 by N matrix.
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
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SROTG
C***FIRST EXECUTABLE STATEMENT  SROTG
      IF (ABS(SA) .LE. ABS(SB)) GO TO 10
C
C *** HERE ABS(SA) .GT. ABS(SB) ***
C
      U = SA + SA
      V = SB / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF SA
C
      R = SQRT(0.25E0 + V**2) * U
C
C     NOTE THAT SC IS POSITIVE
C
      SC = SA / R
      SS = V * (SC + SC)
      SB = SS
      SA = R
      RETURN
C
C *** HERE ABS(SA) .LE. ABS(SB) ***
C
   10 IF (SB .EQ. 0.0E0) GO TO 20
      U = SB + SB
      V = SA / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF SB
C     (R IS IMMEDIATELY STORED IN SA)
C
      SA = SQRT(0.25E0 + V**2) * U
C
C     NOTE THAT SS IS POSITIVE
C
      SS = SB / SA
      SC = V * (SS + SS)
      IF (SC .EQ. 0.0E0) GO TO 15
      SB = 1.0E0 / SC
      RETURN
   15 SB = 1.0E0
      RETURN
C
C *** HERE SA = SB = 0.0 ***
C
   20 SC = 1.0E0
      SS = 0.0E0
      RETURN
C
      END
