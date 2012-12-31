*DECK EXPREL
      FUNCTION EXPREL (X)
C***BEGIN PROLOGUE  EXPREL
C***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      SINGLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
C Taylor series is used.  If X is negative, the reflection formula
C         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
C may be used.  This reflection formula will be of use when the
C evaluation for small ABS(X) is done by Chebyshev series rather than
C Taylor series.  EXPREL and X are single precision.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  EXPREL
      LOGICAL FIRST
      SAVE NTERMS, XBND, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  EXPREL
      IF (FIRST) THEN
         ALNEPS = LOG(R1MACH(3))
         XN = 3.72 - 0.3*ALNEPS
         XLN = LOG((XN+1.0)/1.36)
         NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
         XBND = R1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      ABSX = ABS(X)
      IF (ABSX.GT.0.5) EXPREL = (EXP(X) - 1.0) / X
      IF (ABSX.GT.0.5) RETURN
C
      EXPREL = 1.0
      IF (ABSX.LT.XBND) RETURN
C
      EXPREL = 0.0
      DO 20 I=1,NTERMS
        EXPREL = 1.0 + EXPREL*X/(NTERMS+2-I)
 20   CONTINUE
C
      RETURN
      END
