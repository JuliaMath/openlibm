*DECK DLI
      DOUBLE PRECISION FUNCTION DLI (X)
C***BEGIN PROLOGUE  DLI
C***PURPOSE  Compute the logarithmic integral.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C5
C***TYPE      DOUBLE PRECISION (ALI-S, DLI-D)
C***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLI(X) calculates the double precision logarithmic integral
C for double precision argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DEI, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DLI
      DOUBLE PRECISION X, DEI
C***FIRST EXECUTABLE STATEMENT  DLI
      IF (X .LE. 0.D0) CALL XERMSG ('SLATEC', 'DLI',
     +   'LOG INTEGRAL UNDEFINED FOR X LE 0', 1, 2)
      IF (X .EQ. 1.D0) CALL XERMSG ('SLATEC', 'DLI',
     +   'LOG INTEGRAL UNDEFINED FOR X = 0', 2, 2)
C
      DLI = DEI (LOG(X))
C
      RETURN
      END
