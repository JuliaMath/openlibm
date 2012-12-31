*DECK CACOS
      COMPLEX FUNCTION CACOS (Z)
C***BEGIN PROLOGUE  CACOS
C***PURPOSE  Compute the complex arc cosine.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      COMPLEX (CACOS-C)
C***KEYWORDS  ARC COSINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CACOS(Z) calculates the complex trigonometric arc cosine of Z.
C The result is in units of radians, and the real part is in the
C first or second quadrant.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CASIN
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CACOS
      COMPLEX Z, CASIN
      SAVE PI2
      DATA PI2 /1.5707963267 9489661923E0/
C***FIRST EXECUTABLE STATEMENT  CACOS
      CACOS = PI2 - CASIN (Z)
C
      RETURN
      END
