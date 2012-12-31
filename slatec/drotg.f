*DECK DROTG
      SUBROUTINE DROTG (DA, DB, DC, DS)
C***BEGIN PROLOGUE  DROTG
C***PURPOSE  Construct a plane Givens rotation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      DOUBLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
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
C       DA  double precision scalar
C       DB  double precision scalar
C
C     --Output--
C       DA  double precision result R
C       DB  double precision result Z
C       DC  double precision result
C       DS  double precision result
C
C     Construct the Givens transformation
C
C         ( DC  DS )
C     G = (        ) ,    DC**2 + DS**2 = 1 ,
C         (-DS  DC )
C
C     which zeros the second entry of the 2-vector  (DA,DB)**T .
C
C     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in
C     storage.  The value of DB is overwritten by a value Z which
C     allows DC and DS to be recovered by the following algorithm.
C
C           If Z=1  set  DC=0.0  and  DS=1.0
C           If ABS(Z) .LT. 1  set  DC=SQRT(1-Z**2)  and  DS=Z
C           If ABS(Z) .GT. 1  set  DC=1/Z  and  DS=SQRT(1-DC**2)
C
C     Normally, the subprogram DROT(N,DX,INCX,DY,INCY,DC,DS) will
C     next be called to apply the transformation to a 2 by N matrix.
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
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DROTG
      DOUBLE PRECISION  DA, DB, DC, DS, U, V, R
C***FIRST EXECUTABLE STATEMENT  DROTG
      IF (ABS(DA) .LE. ABS(DB)) GO TO 10
C
C *** HERE ABS(DA) .GT. ABS(DB) ***
C
      U = DA + DA
      V = DB / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF DA
C
      R = SQRT(0.25D0 + V**2) * U
C
C     NOTE THAT DC IS POSITIVE
C
      DC = DA / R
      DS = V * (DC + DC)
      DB = DS
      DA = R
      RETURN
C
C *** HERE ABS(DA) .LE. ABS(DB) ***
C
   10 IF (DB .EQ. 0.0D0) GO TO 20
      U = DB + DB
      V = DA / U
C
C     NOTE THAT U AND R HAVE THE SIGN OF DB
C     (R IS IMMEDIATELY STORED IN DA)
C
      DA = SQRT(0.25D0 + V**2) * U
C
C     NOTE THAT DS IS POSITIVE
C
      DS = DB / DA
      DC = V * (DS + DS)
      IF (DC .EQ. 0.0D0) GO TO 15
      DB = 1.0D0 / DC
      RETURN
   15 DB = 1.0D0
      RETURN
C
C *** HERE DA = DB = 0.0 ***
C
   20 DC = 1.0D0
      DS = 0.0D0
      RETURN
C
      END
