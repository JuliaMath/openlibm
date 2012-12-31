*DECK CROTG
      SUBROUTINE CROTG (CA, CB, C, S)
C***BEGIN PROLOGUE  CROTG
C***PURPOSE  Construct a Givens transformation.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1B10
C***TYPE      COMPLEX (SROTG-S, DROTG-D, CROTG-C)
C***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
C             LINEAR ALGEBRA, VECTOR
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C    Complex Givens transformation
C
C    Construct the Givens transformation
C
C             (C    S)
C       G  =  (      ),  C**2 + ABS(S)**2 =1,
C             (-S   C)
C
C    which zeros the second entry of the complex 2-vector (CA,CB)**T
C
C    The quantity CA/ABS(CA)*NORM(CA,CB) overwrites CA in storage.
C
C    Input:
C        CA (Complex)
C        CB (Complex)
C
C    Output:
C        CA (Complex)      CA/ABS(CA)*NORM(CA,CB)
C        C  (Real)
C        S  (Complex)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CROTG
      COMPLEX CA, CB, S
      REAL C
      REAL NORM, SCALE
      COMPLEX ALPHA
C***FIRST EXECUTABLE STATEMENT  CROTG
      IF (ABS(CA) .EQ. 0.0) THEN
        C = 0.0
        S = (1.0,0.0)
        CA = CB
      ELSE
        SCALE = ABS(CA) + ABS(CB)
        NORM = SCALE * SQRT((ABS(CA/SCALE))**2 + (ABS(CB/SCALE))**2)
        ALPHA = CA /ABS(CA)
        C = ABS(CA) / NORM
        S = ALPHA * CONJG(CB) / NORM
        CA = ALPHA * NORM
      ENDIF
      RETURN
      END
