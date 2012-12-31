*DECK DCHFEV
      SUBROUTINE DCHFEV (X1, X2, F1, F2, D1, D2, NE, XE, FE, NEXT, IERR)
C***BEGIN PROLOGUE  DCHFEV
C***PURPOSE  Evaluate a cubic polynomial given in Hermite form at an
C            array of points.  While designed for use by DPCHFE, it may
C            be useful directly as an evaluator for a piecewise cubic
C            Hermite function in applications, such as graphing, where
C            the interval is known in advance.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3
C***TYPE      DOUBLE PRECISION (CHFEV-S, DCHFEV-D)
C***KEYWORDS  CUBIC HERMITE EVALUATION, CUBIC POLYNOMIAL EVALUATION,
C             PCHIP
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C          DCHFEV:  Cubic Hermite Function EValuator
C
C     Evaluates the cubic polynomial determined by function values
C     F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
C     XE(J), J=1(1)NE.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  NE, NEXT(2), IERR
C        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
C
C        CALL  DCHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
C
C   Parameters:
C
C     X1,X2 -- (input) endpoints of interval of definition of cubic.
C           (Error return if  X1.EQ.X2 .)
C
C     F1,F2 -- (input) values of function at X1 and X2, respectively.
C
C     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real*8 array of points at which the function is to
C           be evaluated.  If any of the XE are outside the interval
C           [X1,X2], a warning error is returned in NEXT.
C
C     FE -- (output) real*8 array of values of the cubic function
C           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
C
C     NEXT -- (output) integer array indicating number of extrapolation
C           points:
C            NEXT(1) = number of evaluation points to left of interval.
C            NEXT(2) = number of evaluation points to right of interval.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if NE.LT.1 .
C              IERR = -2  if X1.EQ.X2 .
C                (The FE-array has not been changed in either case.)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811019  DATE WRITTEN
C   820803  Minor cosmetic changes for release 1.
C   870813  Corrected XERROR calls for d.p. names(s).
C   890206  Corrected XERROR calls.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890703  Corrected category record.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DCHFEV
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DCHFEV to CHFEV wherever it occurs,
C        b. Change the double precision declaration to real, and
C        c. Change the constant ZERO to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  NE, NEXT(2), IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(*), FE(*)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I
      DOUBLE PRECISION  C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA,
     * ZERO
      SAVE ZERO
      DATA  ZERO /0.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DCHFEV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
C
C  INITIALIZE.
C
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = MIN(ZERO, H)
      XMA = MAX(ZERO, H)
C
C  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
C
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
C                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C3 = (DEL1 + DEL2)/H
C                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
C
C  EVALUATION LOOP.
C
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
C          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
C        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'DCHFEV',
     +   'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
      RETURN
C
 5002 CONTINUE
C     X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'DCHFEV', 'INTERVAL ENDPOINTS EQUAL',
     +   IERR, 1)
      RETURN
C------------- LAST LINE OF DCHFEV FOLLOWS -----------------------------
      END
