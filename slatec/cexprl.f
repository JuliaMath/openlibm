*DECK CEXPRL
      COMPLEX FUNCTION CEXPRL (Z)
C***BEGIN PROLOGUE  CEXPRL
C***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      COMPLEX (EXPREL-S, DEXPRL-D, CEXPRL-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  (EXP(Z)-1)/Z .  For small ABS(Z), we use the Taylor
C series.  We could instead use the expression
C        CEXPRL(Z) = (EXP(X)*EXP(I*Y)-1)/Z
C                  = (X*EXPREL(X) * (1 - 2*SIN(Y/2)**2) - 2*SIN(Y/2)**2
C                                    + I*SIN(Y)*(1+X*EXPREL(X))) / Z
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CEXPRL
      COMPLEX Z
      LOGICAL FIRST
      SAVE NTERMS, RBND, FIRST
      DATA FIRST / .TRUE. /
C***FIRST EXECUTABLE STATEMENT  CEXPRL
      IF (FIRST) THEN
         ALNEPS = LOG(R1MACH(3))
         XN = 3.72 - 0.3*ALNEPS
         XLN = LOG((XN+1.0)/1.36)
         NTERMS = XN - (XN*XLN+ALNEPS)/(XLN+1.36) + 1.5
         RBND = R1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      R = ABS(Z)
      IF (R.GT.0.5) CEXPRL = (EXP(Z) - 1.0) / Z
      IF (R.GT.0.5) RETURN
C
      CEXPRL = (1.0, 0.0)
      IF (R.LT.RBND) RETURN
C
      CEXPRL = (0.0, 0.0)
      DO 20 I=1,NTERMS
        CEXPRL = 1.0 + CEXPRL*Z/(NTERMS+2-I)
 20   CONTINUE
C
      RETURN
      END
