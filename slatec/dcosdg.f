*DECK DCOSDG
      DOUBLE PRECISION FUNCTION DCOSDG (X)
C***BEGIN PROLOGUE  DCOSDG
C***PURPOSE  Compute the cosine of an argument in degrees.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      DOUBLE PRECISION (COSDG-S, DCOSDG-D)
C***KEYWORDS  COSINE, DEGREES, ELEMENTARY FUNCTIONS, FNLIB,
C             TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DCOSDG(X) calculates the double precision trigonometric cosine
C for double precision argument X in units of degrees.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DCOSDG
      DOUBLE PRECISION X, RADDEG
      SAVE RADDEG
      DATA RADDEG / 0.0174532925 1994329576 9236907684 886 D0 /
C***FIRST EXECUTABLE STATEMENT  DCOSDG
      DCOSDG = COS (RADDEG*X)
C
      IF (MOD(X,90.D0).NE.0.D0) RETURN
      N = ABS(X)/90.D0 + 0.5D0
      N = MOD (N, 2)
      IF (N.EQ.0) DCOSDG = SIGN (1.0D0, DCOSDG)
      IF (N.EQ.1) DCOSDG = 0.0D0
C
      RETURN
      END
