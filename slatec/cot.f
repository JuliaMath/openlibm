*DECK COT
      FUNCTION COT (X)
C***BEGIN PROLOGUE  COT
C***PURPOSE  Compute the cotangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      SINGLE PRECISION (COT-S, DCOT-D, CCOT-C)
C***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C COT(X) calculates the cotangent of the real argument X.  X is in
C units of radians.
C
C Series for COT        on the interval  0.          to  6.25000D-02
C                                        with weighted error   3.76E-17
C                                         log weighted error  16.42
C                               significant figures required  15.51
C                                    decimal places required  16.88
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  COT
      DIMENSION COTCS(8)
      LOGICAL FIRST
      SAVE COTCS, PI2REC, NTERMS, XMAX, XSML, XMIN, SQEPS, FIRST
      DATA COTCS( 1) /    .2402591609 8295630E0 /
      DATA COTCS( 2) /   -.0165330316 01500228E0 /
      DATA COTCS( 3) /   -.0000429983 91931724E0 /
      DATA COTCS( 4) /   -.0000001592 83223327E0 /
      DATA COTCS( 5) /   -.0000000006 19109313E0 /
      DATA COTCS( 6) /   -.0000000000 02430197E0 /
      DATA COTCS( 7) /   -.0000000000 00009560E0 /
      DATA COTCS( 8) /   -.0000000000 00000037E0 /
      DATA PI2REC / .01161977236 75813430 E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  COT
      IF (FIRST) THEN
         NTERMS = INITS (COTCS, 8, 0.1*R1MACH(3))
         XMAX = 1.0/R1MACH(4)
         XSML = SQRT (3.0*R1MACH(3))
         XMIN = EXP ( MAX(LOG(R1MACH(1)), -LOG(R1MACH(2))) + 0.01)
         SQEPS = SQRT (R1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (ABS(X) .LT. XMIN) CALL XERMSG ('SLATEC', 'COT',
     +   'ABS(X) IS ZERO OR SO SMALL COT OVERFLOWS', 2, 2)
      IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'COT',
     +   'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 3, 2)
C
C CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
C = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
C = AINT(.625*Y) + AINT(Z) + REM(Z)
C
      AINTY = AINT (Y)
      YREM = Y - AINTY
      PRODBG = 0.625*AINTY
      AINTY = AINT (PRODBG)
      Y = (PRODBG-AINTY) + 0.625*YREM + Y*PI2REC
      AINTY2 = AINT (Y)
      AINTY = AINTY + AINTY2
      Y = Y - AINTY2
C
      IFN = MOD (AINTY, 2.)
      IF (IFN.EQ.1) Y = 1.0 - Y
C
      IF (ABS(X) .GT. 0.5 .AND. Y .LT. ABS(X)*SQEPS) CALL XERMSG
     +   ('SLATEC', 'COT',
     +   'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI ' //
     +   '(N.NE.0)' , 1, 1)
C
      IF (Y.GT.0.25) GO TO 20
      COT = 1.0/X
      IF (Y.GT.XSML) COT = (0.5 + CSEVL (32.0*Y*Y-1., COTCS, NTERMS)) /Y
      GO TO 40
C
 20   IF (Y.GT.0.5) GO TO 30
      COT = (0.5 + CSEVL (8.0*Y*Y-1., COTCS, NTERMS)) / (0.5*Y)
      COT = (COT**2 - 1.0) * 0.5 / COT
      GO TO 40
C
 30   COT = (0.5 + CSEVL (2.0*Y*Y-1., COTCS, NTERMS)) / (0.25*Y)
      COT = (COT**2 - 1.0) * 0.5 / COT
      COT = (COT**2 - 1.0) * 0.5 / COT
C
 40   IF (X.NE.0.) COT = SIGN (COT, X)
      IF (IFN.EQ.1) COT = -COT
C
      RETURN
      END
