*DECK CCBRT
      COMPLEX FUNCTION CCBRT (Z)
C***BEGIN PROLOGUE  CCBRT
C***PURPOSE  Compute the cube root.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C2
C***TYPE      COMPLEX (CBRT-S, DCBRT-D, CCBRT-C)
C***KEYWORDS  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CCBRT(Z) calculates the complex cube root of Z.  The principal root
C for which -PI .LT. arg(Z) .LE. +PI is returned.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CARG, CBRT
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CCBRT
      COMPLEX Z
C***FIRST EXECUTABLE STATEMENT  CCBRT
      THETA = CARG(Z) / 3.0
      R = CBRT (ABS(Z))
C
      CCBRT = CMPLX (R*COS(THETA), R*SIN(THETA))
C
      RETURN
      END
