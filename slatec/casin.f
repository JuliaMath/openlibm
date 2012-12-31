*DECK CASIN
      COMPLEX FUNCTION CASIN (ZINP)
C***BEGIN PROLOGUE  CASIN
C***PURPOSE  Compute the complex arc sine.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      COMPLEX (CASIN-C)
C***KEYWORDS  ARC SINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CASIN(ZINP) calculates the complex trigonometric arc sine of ZINP.
C The result is in units of radians, and the real part is in the first
C or fourth quadrant.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  CASIN
      COMPLEX ZINP, Z, Z2, SQZP1, CI
      LOGICAL FIRST
      SAVE PI2, PI, CI, NTERMS, RMIN, FIRST
      DATA PI2 /1.5707963267 9489661923E0/
      DATA PI /3.1415926535 8979324E0/
      DATA CI /(0.,1.)/
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  CASIN
      IF (FIRST) THEN
C NTERMS = LOG(EPS)/LOG(RMAX)  WHERE RMAX = 0.1
         NTERMS = -0.4343*LOG(R1MACH(3))
         RMIN = SQRT (6.0*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Z = ZINP
      R = ABS (Z)
      IF (R.GT.0.1) GO TO 30
C
      CASIN = Z
      IF (R.LT.RMIN) RETURN
C
      CASIN = (0.0, 0.0)
      Z2 = Z*Z
      DO 20 I=1,NTERMS
        TWOI = 2*(NTERMS-I) + 1
        CASIN = 1.0/TWOI + TWOI*CASIN*Z2/(TWOI+1.0)
 20   CONTINUE
      CASIN = Z*CASIN
      RETURN
C
 30   IF (REAL(ZINP).LT.0.0) Z = -ZINP
C
      SQZP1 = SQRT (Z+1.0)
      IF (AIMAG(SQZP1).LT.0.) SQZP1 = -SQZP1
      CASIN = PI2 - CI * LOG (Z + SQZP1*SQRT(Z-1.0))
C
      IF (REAL(CASIN).GT.PI2) CASIN = PI - CASIN
      IF (REAL(CASIN).LE.(-PI2)) CASIN = -PI - CASIN
      IF (REAL(ZINP).LT.0.) CASIN = -CASIN
C
      RETURN
      END
