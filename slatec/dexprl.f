*DECK DEXPRL
      DOUBLE PRECISION FUNCTION DEXPRL (X)
C***BEGIN PROLOGUE  DEXPRL
C***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      DOUBLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the
C Taylor series is used.  If X is negative the reflection formula
C         EXPREL(X) = EXP(X) * EXPREL(ABS(X))
C may be used.  This reflection formula will be of use when the
C evaluation for small ABS(X) is done by Chebyshev series rather than
C Taylor series.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DEXPRL
      DOUBLE PRECISION X, ABSX, ALNEPS, XBND, XLN, XN,  D1MACH
      LOGICAL FIRST
      SAVE NTERMS, XBND, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DEXPRL
      IF (FIRST) THEN
         ALNEPS = LOG(D1MACH(3))
         XN = 3.72D0 - 0.3D0*ALNEPS
         XLN = LOG((XN+1.0D0)/1.36D0)
         NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36D0) + 1.5D0
         XBND = D1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      ABSX = ABS(X)
      IF (ABSX.GT.0.5D0) DEXPRL = (EXP(X)-1.0D0)/X
      IF (ABSX.GT.0.5D0) RETURN
C
      DEXPRL = 1.0D0
      IF (ABSX.LT.XBND) RETURN
C
      DEXPRL = 0.0D0
      DO 20 I=1,NTERMS
        DEXPRL = 1.0D0 + DEXPRL*X/(NTERMS+2-I)
 20   CONTINUE
C
      RETURN
      END
