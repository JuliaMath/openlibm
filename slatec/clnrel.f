*DECK CLNREL
      COMPLEX FUNCTION CLNREL (Z)
C***BEGIN PROLOGUE  CLNREL
C***PURPOSE  Evaluate ln(1+X) accurate in the sense of relative error.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      COMPLEX (ALNREL-S, DLNREL-D, CLNREL-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CLNREL(Z) = LOG(1+Z) with relative error accuracy near Z = 0.
C Let   RHO = ABS(Z)  and
C       R**2 = ABS(1+Z)**2 = (1+X)**2 + Y**2 = 1 + 2*X + RHO**2 .
C Now if RHO is small we may evaluate CLNREL(Z) accurately by
C       LOG(1+Z) = CMPLX  (LOG(R), CARG(1+Z))
C                 = CMPLX  (0.5*LOG(R**2), CARG(1+Z))
C                 = CMPLX  (0.5*ALNREL(2*X+RHO**2), CARG(1+Z))
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ALNREL, CARG, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  CLNREL
      COMPLEX Z
      SAVE SQEPS
      DATA SQEPS /0.0/
C***FIRST EXECUTABLE STATEMENT  CLNREL
      IF (SQEPS.EQ.0.) SQEPS = SQRT (R1MACH(4))
C
      IF (ABS(1.+Z) .LT. SQEPS) CALL XERMSG ('SLATEC', 'CLNREL',
     +   'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR -1', 1, 1)
C
      RHO = ABS(Z)
      IF (RHO.GT.0.375) CLNREL = LOG (1.0+Z)
      IF (RHO.GT.0.375) RETURN
C
      X = REAL(Z)
      CLNREL = CMPLX (0.5*ALNREL(2.*X+RHO**2), CARG(1.0+Z))
C
      RETURN
      END
