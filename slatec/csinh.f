*DECK CSINH
      COMPLEX FUNCTION CSINH (Z)
C***BEGIN PROLOGUE  CSINH
C***PURPOSE  Compute the complex hyperbolic sine.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4C
C***TYPE      COMPLEX (CSINH-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC SINE
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CSINH(Z) calculates the complex hyperbolic sine of complex
C argument Z.  Z is in units of radians.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CSINH
      COMPLEX Z, CI
      SAVE CI
      DATA CI /(0.,1.)/
C***FIRST EXECUTABLE STATEMENT  CSINH
      CSINH = -CI*SIN(CI*Z)
C
      RETURN
      END
