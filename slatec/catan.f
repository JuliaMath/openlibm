*DECK CATAN
      COMPLEX FUNCTION CATAN (Z)
C***BEGIN PROLOGUE  CATAN
C***PURPOSE  Compute the complex arc tangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      COMPLEX (CATAN-C)
C***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C CATAN(Z) calculates the complex trigonometric arc tangent of Z.
C The result is in units of radians, and the real part is in the first
C or fourth quadrant.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770801  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  CATAN
      COMPLEX Z, Z2
      LOGICAL FIRST
      SAVE PI2, NTERMS, SQEPS, RMIN, RMAX, FIRST
      DATA PI2 / 1.5707963267 9489661923E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  CATAN
      IF (FIRST) THEN
C NTERMS = LOG(EPS)/LOG(RBND) WHERE RBND = 0.1
         NTERMS = -0.4343*LOG(R1MACH(3)) + 1.0
         SQEPS = SQRT(R1MACH(4))
         RMIN = SQRT (3.0*R1MACH(3))
         RMAX = 1.0/R1MACH(3)
      ENDIF
      FIRST = .FALSE.
C
      R = ABS(Z)
      IF (R.GT.0.1) GO TO 30
C
      CATAN = Z
      IF (R.LT.RMIN) RETURN
C
      CATAN = (0.0, 0.0)
      Z2 = Z*Z
      DO 20 I=1,NTERMS
        TWOI = 2*(NTERMS-I) + 1
        CATAN = 1.0/TWOI - Z2*CATAN
 20   CONTINUE
      CATAN = Z*CATAN
      RETURN
C
 30   IF (R.GT.RMAX) GO TO 50
      X = REAL(Z)
      Y = AIMAG(Z)
      R2 = R*R
      IF (R2 .EQ. 1.0 .AND. X .EQ. 0.0) CALL XERMSG ('SLATEC', 'CATAN',
     +   'Z IS +I OR -I', 2, 2)
      IF (ABS(R2-1.0).GT.SQEPS) GO TO 40
      IF (ABS(CMPLX(1.0, 0.0)+Z*Z) .LT. SQEPS) CALL XERMSG ('SLATEC',
     +   'CATAN', 'ANSWER LT HALF PRECISION, Z**2 CLOSE TO -1', 1, 1)
C
 40   XANS = 0.5*ATAN2(2.0*X, 1.0-R2)
      YANS = 0.25*LOG((R2+2.0*Y+1.0)/(R2-2.0*Y+1.0))
      CATAN = CMPLX (XANS, YANS)
      RETURN
C
 50   CATAN = CMPLX (PI2, 0.)
      IF (REAL(Z).LT.0.0) CATAN = CMPLX(-PI2,0.0)
      RETURN
C
      END
