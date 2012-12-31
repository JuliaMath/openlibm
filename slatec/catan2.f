*DECK CATAN2
      COMPLEX FUNCTION CATAN2 (CSN, CCS)
C***BEGIN PROLOGUE  CATAN2
C***PURPOSE  Compute the complex arc tangent in the proper quadrant.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      COMPLEX (CATAN2-C)
C***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, POLAR ANGEL,
C             QUADRANT, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CATAN2(CSN,CCS) calculates the complex trigonometric arc
C tangent of the ratio CSN/CCS and returns a result whose real
C part is in the correct quadrant (within a multiple of 2*PI).  The
C result is in units of radians and the real part is between -PI
C and +PI.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CATAN, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  CATAN2
      COMPLEX CSN, CCS, CATAN
      SAVE PI
      DATA PI / 3.1415926535 8979323846E0 /
C***FIRST EXECUTABLE STATEMENT  CATAN2
      IF (ABS(CCS).EQ.0.) GO TO 10
C
      CATAN2 = CATAN (CSN/CCS)
      IF (REAL(CCS).LT.0.) CATAN2 = CATAN2 + PI
      IF (REAL(CATAN2).GT.PI) CATAN2 = CATAN2 - 2.0*PI
      RETURN
C
 10   IF (ABS(CSN) .EQ. 0.) CALL XERMSG ('SLATEC', 'CATAN2',
     +   'CALLED WITH BOTH ARGUMENTS ZERO', 1, 2)
C
      CATAN2 = CMPLX (SIGN(0.5*PI,REAL(CSN)), 0.0)
C
      RETURN
      END
