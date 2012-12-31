*DECK C9LN2R
      COMPLEX FUNCTION C9LN2R (Z)
C***BEGIN PROLOGUE  C9LN2R
C***SUBSIDIARY
C***PURPOSE  Evaluate LOG(1+Z) from second order relative accuracy so
C            that  LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      COMPLEX (R9LN2R-S, D9LN2R-D, C9LN2R-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  LOG(1+Z)  from 2-nd order with relative error accuracy so
C that     LOG(1+Z) = Z - Z**2/2 + Z**3*C9LN2R(Z).
C
C Now  LOG(1+Z) = 0.5*LOG(1+2*X+ABS(Z)**2) + I*CARG(1+Z),
C where X = REAL(Z)  and  Y = AIMAG(Z).
C We find
C     Z**3 * C9LN2R(Z) = -X*ABS(Z)**2 - 0.25*ABS(Z)**4
C        + (2*X+ABS(Z)**2)**3 * R9LN2R(2*X+ABS(Z)**2)
C        + I * (CARG(1+Z) + (X-1)*Y)
C The imaginary part must be evaluated carefully as
C     (ATAN(Y/(1+X)) - Y/(1+X)) + Y/(1+X) - (1-X)*Y
C       = (Y/(1+X))**3 * R9ATN1(Y/(1+X)) + X**2*Y/(1+X)
C
C Now we divide through by Z**3 carefully.  Write
C     1/Z**3 = (X-I*Y)/ABS(Z)**3 * (1/ABS(Z)**3)
C then   C9LN2R(Z) = ((X-I*Y)/ABS(Z))**3 * (-X/ABS(Z) - ABS(Z)/4
C        + 0.5*((2*X+ABS(Z)**2)/ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2)
C        + I*Y/(ABS(Z)*(1+X)) * ((X/ABS(Z))**2 +
C          + (Y/(ABS(Z)*(1+X)))**2 * R9ATN1(Y/(1+X)) ) )
C
C If we let  XZ = X/ABS(Z)  and  YZ = Y/ABS(Z)  we may write
C     C9LN2R(Z) = (XZ-I*YZ)**3 * (-XZ - ABS(Z)/4
C        + 0.5*(2*XZ+ABS(Z))**3 * R9LN2R(2*X+ABS(Z)**2)
C        + I*YZ/(1+X) * (XZ**2 + (YZ/(1+X))**2*R9ATN1(Y/(1+X)) ))
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R9ATN1, R9LN2R
C***REVISION HISTORY  (YYMMDD)
C   780401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  C9LN2R
      COMPLEX Z
C***FIRST EXECUTABLE STATEMENT  C9LN2R
      X = REAL (Z)
      Y = AIMAG (Z)
C
      CABSZ = ABS(Z)
      IF (CABSZ.GT.0.8125) GO TO 20
C
      C9LN2R = CMPLX (1.0/3.0, 0.0)
      IF (CABSZ.EQ.0.0) RETURN
C
      XZ = X/CABSZ
      YZ = Y/CABSZ
C
      ARG = 2.0*XZ + CABSZ
      RPART = 0.5*ARG**3*R9LN2R(CABSZ*ARG) - XZ - 0.25*CABSZ
      Y1X = YZ/(1.0+X)
      AIPART = Y1X * (XZ**2 + Y1X**2*R9ATN1(CABSZ*Y1X) )
C
      C9LN2R = CMPLX(XZ,-YZ)**3 * CMPLX(RPART,AIPART)
      RETURN
C
 20   C9LN2R = (LOG(1.0+Z) - Z*(1.0-0.5*Z)) / Z**3
      RETURN
C
      END
