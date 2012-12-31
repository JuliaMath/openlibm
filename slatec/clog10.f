*DECK CLOG10
      COMPLEX FUNCTION CLOG10 (Z)
C***BEGIN PROLOGUE  CLOG10
C***PURPOSE  Compute the principal value of the complex base 10
C            logarithm.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      COMPLEX (CLOG10-C)
C***KEYWORDS  BASE TEN LOGARITHM, ELEMENTARY FUNCTIONS, FNLIB
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CLOG10(Z) calculates the principal value of the complex common
C or base 10 logarithm of Z for -PI .LT. arg(Z) .LE. +PI.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CLOG10
      COMPLEX Z
      SAVE ALOGE
      DATA ALOGE / 0.4342944819 0325182765E0 /
C***FIRST EXECUTABLE STATEMENT  CLOG10
      CLOG10 = ALOGE * LOG(Z)
C
      RETURN
      END
