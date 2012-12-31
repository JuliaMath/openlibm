*DECK DPPVAL
      DOUBLE PRECISION FUNCTION DPPVAL (LDC, C, XI, LXI, K, IDERIV, X,
     +   INPPV)
C***BEGIN PROLOGUE  DPPVAL
C***PURPOSE  Calculate the value of the IDERIV-th derivative of the
C            B-spline from the PP-representation.
C***LIBRARY   SLATEC
C***CATEGORY  E3, K6
C***TYPE      DOUBLE PRECISION (PPVAL-S, DPPVAL-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract    **** a double precision routine ****
C         DPPVAL is the PPVALU function of the reference.
C
C         DPPVAL calculates (at X) the value of the IDERIV-th
C         derivative of the B-spline from the PP-representation
C         (C,XI,LXI,K).  The Taylor expansion about XI(J) for X in
C         the interval XI(J) .LE. X .LT. XI(J+1) is evaluated, J=1,LXI.
C         Right limiting values at X=XI(J) are obtained.  DPPVAL will
C         extrapolate beyond XI(1) and XI(LXI+1).
C
C         To obtain left limiting values (left derivatives) at XI(J)
C         replace LXI by J-1 and set X=XI(J),J=2,LXI+1.
C
C     Description of Arguments
C
C         Input      C,XI,X are double precision
C          LDC     - leading dimension of C matrix, LDC .GE. K
C          C       - matrix of dimension at least (K,LXI) containing
C                    right derivatives at break points XI(*).
C          XI      - break point vector of length LXI+1
C          LXI     - number of polynomial pieces
C          K       - order of B-spline, K .GE. 1
C          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
C                    IDERIV=0 gives the B-spline value
C          X       - argument, XI(1) .LE. X .LE. XI(LXI+1)
C          INPPV   - an initialization parameter which must be set
C                    to 1 the first time DPPVAL is called.
C
C         Output     DPPVAL is double precision
C          INPPV   - INPPV contains information for efficient process-
C                    ing after the initial call and INPPV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INPPV parameters.
C          DPPVAL  - value of the IDERIV-th derivative at X
C
C     Error Conditions
C         Improper input is a fatal error
C
C***REFERENCES  Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C***ROUTINES CALLED  DINTRV, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPPVAL
C
      INTEGER I, IDERIV, INPPV, J, K, LDC, LXI, NDUMMY, KK
      DOUBLE PRECISION C, DX, X, XI
      DIMENSION XI(*), C(LDC,*)
C***FIRST EXECUTABLE STATEMENT  DPPVAL
      DPPVAL = 0.0D0
      IF(K.LT.1) GO TO 90
      IF(LDC.LT.K) GO TO 80
      IF(LXI.LT.1) GO TO 85
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 95
      I = K - IDERIV
      KK = I
      CALL DINTRV(XI, LXI, X, INPPV, I, NDUMMY)
      DX = X - XI(I)
      J = K
   10 DPPVAL = (DPPVAL/KK)*DX + C(J,I)
      J = J - 1
      KK = KK - 1
      IF (KK.GT.0) GO TO 10
      RETURN
C
C
   80 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'LDC DOES NOT SATISFY LDC.GE.K',
     +   2, 1)
      RETURN
   85 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'LXI DOES NOT SATISFY LXI.GE.1',
     +   2, 1)
      RETURN
   90 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
   95 CONTINUE
      CALL XERMSG ('SLATEC', 'DPPVAL',
     +   'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K', 2, 1)
      RETURN
      END
