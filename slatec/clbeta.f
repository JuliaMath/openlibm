*DECK CLBETA
      COMPLEX FUNCTION CLBETA (A, B)
C***BEGIN PROLOGUE  CLBETA
C***PURPOSE  Compute the natural logarithm of the complete Beta
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7B
C***TYPE      COMPLEX (ALBETA-S, DLBETA-D, CLBETA-C)
C***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CLBETA computes the natural log of the complex valued complete beta
C function of complex parameters A and B.  This is a preliminary version
C which is not accurate.
C
C Input Parameters:
C       A   complex and the real part of A positive
C       B   complex and the real part of B positive
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CLNGAM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  CLBETA
      COMPLEX A, B, CLNGAM
C***FIRST EXECUTABLE STATEMENT  CLBETA
      IF (REAL(A) .LE. 0.0 .OR. REAL(B) .LE. 0.0) CALL XERMSG ('SLATEC',
     +   'CLBETA', 'REAL PART OF BOTH ARGUMENTS MUST BE GT 0', 1, 2)
C
      CLBETA = CLNGAM(A) + CLNGAM(B) - CLNGAM(A+B)
C
      RETURN
      END
