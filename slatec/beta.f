*DECK BETA
      FUNCTION BETA (A, B)
C***BEGIN PROLOGUE  BETA
C***PURPOSE  Compute the complete Beta function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7B
C***TYPE      SINGLE PRECISION (BETA-S, DBETA-D, CBETA-C)
C***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BETA computes the complete beta function.
C
C Input Parameters:
C       A   real and positive
C       B   real and positive
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALBETA, GAMLIM, GAMMA, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  BETA
      EXTERNAL GAMMA
      SAVE XMAX, ALNSML
      DATA XMAX, ALNSML /0., 0./
C***FIRST EXECUTABLE STATEMENT  BETA
      IF (ALNSML.NE.0.0) GO TO 10
      CALL GAMLIM (XMIN, XMAX)
      ALNSML = LOG(R1MACH(1))
C
 10   IF (A .LE. 0. .OR. B .LE. 0.) CALL XERMSG ('SLATEC', 'BETA',
     +   'BOTH ARGUMENTS MUST BE GT 0', 2, 2)
C
      IF (A+B.LT.XMAX) BETA = GAMMA(A) * GAMMA(B) / GAMMA(A+B)
      IF (A+B.LT.XMAX) RETURN
C
      BETA = ALBETA (A, B)
      IF (BETA .LT. ALNSML) CALL XERMSG ('SLATEC', 'BETA',
     +   'A AND/OR B SO BIG BETA UNDERFLOWS', 1, 2)
C
      BETA = EXP (BETA)
C
      RETURN
      END
