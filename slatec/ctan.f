*DECK CTAN
      COMPLEX FUNCTION CTAN (Z)
C***BEGIN PROLOGUE  CTAN
C***PURPOSE  Compute the complex tangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      COMPLEX (CTAN-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, TANGENT, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CTAN(Z) calculates the complex trigonometric tangent of complex
C argument Z.  Z is in units of radians.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERCLR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  CTAN
      COMPLEX Z
      SAVE SQEPS
      DATA SQEPS /0./
C***FIRST EXECUTABLE STATEMENT  CTAN
      IF (SQEPS.EQ.0.) SQEPS = SQRT (R1MACH(4))
C
      X2 = 2.0*REAL(Z)
      Y2 = 2.0*AIMAG(Z)
C
      SN2X = SIN (X2)
      CALL XERCLR
C
      DEN = COS(X2) + COSH(Y2)
      IF (DEN .EQ. 0.) CALL XERMSG ('SLATEC', 'CTAN',
     +   'TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR 3*PI/2 AND Y IS 0)',
     +   2, 2)
C
      IF (ABS(DEN).GT.MAX(ABS(X2),1.)*SQEPS) GO TO 10
      CALL XERCLR
      CALL XERMSG ('SLATEC', 'CTAN',
     +   'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR ' //
     +   'PI/2 OR 3*PI/2', 1, 1)
C
 10   CTAN = CMPLX (SN2X/DEN, SINH(Y2)/DEN)
C
      RETURN
      END
