*DECK C0LGMC
      COMPLEX FUNCTION C0LGMC (Z)
C***BEGIN PROLOGUE  C0LGMC
C***PURPOSE  Evaluate (Z+0.5)*LOG((Z+1.)/Z) - 1.0 with relative
C            accuracy.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      COMPLEX (C0LGMC-C)
C***KEYWORDS  FNLIB, GAMMA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  (Z+0.5)*LOG((Z+1.0)/Z) - 1.0  with relative error accuracy
C Let Q = 1.0/Z so that
C     (Z+0.5)*LOG(1+1/Z) - 1 = (Z+0.5)*(LOG(1+Q) - Q + Q*Q/2) - Q*Q/4
C        = (Z+0.5)*Q**3*C9LN2R(Q) - Q**2/4,
C where  C9LN2R  is (LOG(1+Q) - Q + 0.5*Q**2) / Q**3.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  C9LN2R, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   780401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  C0LGMC
      COMPLEX Z, Q, C9LN2R
      SAVE RBIG
      DATA RBIG / 0.0 /
C***FIRST EXECUTABLE STATEMENT  C0LGMC
      IF (RBIG.EQ.0.0) RBIG = 1.0/R1MACH(3)
C
      CABSZ = ABS(Z)
      IF (CABSZ.GT.RBIG) C0LGMC = -(Z+0.5)*LOG(Z) - Z
      IF (CABSZ.GT.RBIG) RETURN
C
      Q = 1.0/Z
      IF (CABSZ.LE.1.23) C0LGMC = (Z+0.5)*LOG(1.0+Q) - 1.0
      IF (CABSZ.GT.1.23) C0LGMC = ((1.+.5*Q)*C9LN2R(Q) - .25) * Q**2
C
      RETURN
      END
