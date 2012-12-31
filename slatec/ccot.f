*DECK CCOT
      COMPLEX FUNCTION CCOT (Z)
C***BEGIN PROLOGUE  CCOT
C***PURPOSE  Compute the cotangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      COMPLEX (COT-S, DCOT-D, CCOT-C)
C***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CCOT(Z) calculates the complex trigonometric cotangent of Z.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERCLR, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  CCOT
      COMPLEX Z
      SAVE SQEPS
      DATA SQEPS /0./
C***FIRST EXECUTABLE STATEMENT  CCOT
      IF (SQEPS.EQ.0.) SQEPS = SQRT (R1MACH(4))
C
      X2 = 2.0*REAL(Z)
      Y2 = 2.0*AIMAG(Z)
C
      SN2X = SIN (X2)
      CALL XERCLR
C
      DEN = COSH(Y2) - COS(X2)
      IF (DEN .EQ. 0.) CALL XERMSG ('SLATEC', 'CCOT',
     +   'COT IS SINGULAR FOR INPUT Z (X IS 0 OR PI AND Y IS 0)', 2, 2)
C
      IF (ABS(DEN).GT.MAX(ABS(X2),1.)*SQEPS) GO TO 10
      CALL XERCLR
      CALL XERMSG ('SLATEC', 'CCOT',
     +   'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR ' //
     +   '0 OR PI', 1, 1)
C
 10   CCOT = CMPLX (SN2X/DEN, -SINH(Y2)/DEN)
C
      RETURN
      END
