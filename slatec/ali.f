*DECK ALI
      FUNCTION ALI (X)
C***BEGIN PROLOGUE  ALI
C***PURPOSE  Compute the logarithmic integral.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C5
C***TYPE      SINGLE PRECISION (ALI-S, DLI-D)
C***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C ALI(X) computes the logarithmic integral; i.e., the
C integral from 0.0 to X of (1.0/ln(t))dt.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  EI, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  ALI
C***FIRST EXECUTABLE STATEMENT  ALI
      IF (X .LE. 0.0) CALL XERMSG ('SLATEC', 'ALI',
     +   'LOG INTEGRAL UNDEFINED FOR X LE 0', 1, 2)
      IF (X .EQ. 1.0) CALL XERMSG ('SLATEC', 'ALI',
     +   'LOG INTEGRAL UNDEFINED FOR X = 1', 2, 2)
C
      ALI = EI (LOG(X) )
C
      RETURN
      END
