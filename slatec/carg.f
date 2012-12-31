*DECK CARG
      FUNCTION CARG (Z)
C***BEGIN PROLOGUE  CARG
C***PURPOSE  Compute the argument of a complex number.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  A4A
C***TYPE      COMPLEX (CARG-C)
C***KEYWORDS  ARGUMENT OF A COMPLEX NUMBER, ELEMENTARY FUNCTIONS, FNLIB
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CARG(Z) calculates the argument of the complex number Z.  Note
C that CARG returns a real result.  If Z = X+iY, then CARG is ATAN(Y/X),
C except when both X and Y are zero, in which case the result
C will be zero.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CARG
      COMPLEX Z
C***FIRST EXECUTABLE STATEMENT  CARG
      CARG = 0.0
      IF (REAL(Z).NE.0. .OR. AIMAG(Z).NE.0.) CARG =
     1  ATAN2 (AIMAG(Z), REAL(Z))
C
      RETURN
      END
