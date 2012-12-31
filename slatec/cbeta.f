*DECK CBETA
      COMPLEX FUNCTION CBETA (A, B)
C***BEGIN PROLOGUE  CBETA
C***PURPOSE  Compute the complete Beta function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7B
C***TYPE      COMPLEX (BETA-S, DBETA-D, CBETA-C)
C***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CBETA computes the complete beta function of complex parameters A
C and B.
C Input Parameters:
C       A   complex and the real part of A positive
C       B   complex and the real part of B positive
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CGAMMA, CLBETA, GAMLIM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  CBETA
      COMPLEX A, B, CGAMMA, CLBETA
      EXTERNAL CGAMMA
      SAVE XMAX
      DATA XMAX / 0.0 /
C***FIRST EXECUTABLE STATEMENT  CBETA
      IF (XMAX.EQ.0.0) THEN
         CALL GAMLIM (XMIN, XMAXT)
         XMAX = XMAXT
      ENDIF
C
      IF (REAL(A) .LE. 0.0 .OR. REAL(B) .LE. 0.0) CALL XERMSG ('SLATEC',
     +   'CBETA', 'REAL PART OF BOTH ARGUMENTS MUST BE GT 0', 1, 2)
C
      IF (REAL(A)+REAL(B).LT.XMAX) CBETA = CGAMMA(A) * (CGAMMA(B)/
     1  CGAMMA(A+B) )
      IF (REAL(A)+REAL(B).LT.XMAX) RETURN
C
      CBETA = EXP (CLBETA(A, B))
C
      RETURN
      END
