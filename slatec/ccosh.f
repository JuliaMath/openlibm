*DECK CCOSH
      COMPLEX FUNCTION CCOSH (Z)
C***BEGIN PROLOGUE  CCOSH
C***PURPOSE  Compute the complex hyperbolic cosine.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4C
C***TYPE      COMPLEX (CCOSH-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC COSINE
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CCOSH(Z) calculates the complex hyperbolic cosine of Z.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CCOSH
      COMPLEX Z, CI
      SAVE CI
      DATA CI /(0.,1.)/
C***FIRST EXECUTABLE STATEMENT  CCOSH
      CCOSH = COS (CI*Z)
C
      RETURN
      END
