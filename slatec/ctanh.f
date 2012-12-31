*DECK CTANH
      COMPLEX FUNCTION CTANH (Z)
C***BEGIN PROLOGUE  CTANH
C***PURPOSE  Compute the complex hyperbolic tangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4C
C***TYPE      COMPLEX (CTANH-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC TANGENT
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CTANH(Z) calculates the complex hyperbolic tangent of complex
C argument Z.  Z is in units of radians.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CTAN
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CTANH
      COMPLEX Z, CI, CTAN
      SAVE CI
      DATA CI /(0.,1.)/
C***FIRST EXECUTABLE STATEMENT  CTANH
      CTANH = -CI*CTAN(CI*Z)
C
      RETURN
      END
