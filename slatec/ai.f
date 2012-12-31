*DECK AI
      FUNCTION AI (X)
C***BEGIN PROLOGUE  AI
C***PURPOSE  Evaluate the Airy function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10D
C***TYPE      SINGLE PRECISION (AI-S, DAI-D)
C***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C AI(X) computes the Airy function Ai(X)
C Series for AIF        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.09E-19
C                                         log weighted error  18.96
C                               significant figures required  17.76
C                                    decimal places required  19.44
C
C Series for AIG        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   1.51E-17
C                                         log weighted error  16.82
C                               significant figures required  15.19
C                                    decimal places required  17.27
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  AIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  AI
      DIMENSION AIFCS(9), AIGCS(8)
      LOGICAL FIRST
      SAVE AIFCS, AIGCS, NAIF, NAIG, X3SML, XMAX, FIRST
      DATA AIFCS( 1) /   -.0379713584 9666999750E0 /
      DATA AIFCS( 2) /    .0591918885 3726363857E0 /
      DATA AIFCS( 3) /    .0009862928 0577279975E0 /
      DATA AIFCS( 4) /    .0000068488 4381907656E0 /
      DATA AIFCS( 5) /    .0000000259 4202596219E0 /
      DATA AIFCS( 6) /    .0000000000 6176612774E0 /
      DATA AIFCS( 7) /    .0000000000 0010092454E0 /
      DATA AIFCS( 8) /    .0000000000 0000012014E0 /
      DATA AIFCS( 9) /    .0000000000 0000000010E0 /
      DATA AIGCS( 1) /    .0181523655 8116127E0 /
      DATA AIGCS( 2) /    .0215725631 6601076E0 /
      DATA AIGCS( 3) /    .0002567835 6987483E0 /
      DATA AIGCS( 4) /    .0000014265 2141197E0 /
      DATA AIGCS( 5) /    .0000000045 7211492E0 /
      DATA AIGCS( 6) /    .0000000000 0952517E0 /
      DATA AIGCS( 7) /    .0000000000 0001392E0 /
      DATA AIGCS( 8) /    .0000000000 0000001E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  AI
      IF (FIRST) THEN
         NAIF = INITS (AIFCS, 9, 0.1*R1MACH(3))
         NAIG = INITS (AIGCS, 8, 0.1*R1MACH(3))
C
         X3SML = R1MACH(3)**0.3334
         XMAXT = (-1.5*LOG(R1MACH(1)))**0.6667
         XMAX = XMAXT - XMAXT*LOG(XMAXT)/
     *                   (4.0*SQRT(XMAXT)+1.0) - 0.01
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GE.(-1.0)) GO TO 20
      CALL R9AIMP (X, XM, THETA)
      AI = XM * COS(THETA)
      RETURN
C
 20   IF (X.GT.1.0) GO TO 30
      Z = 0.0
      IF (ABS(X).GT.X3SML) Z = X**3
      AI = 0.375 + (CSEVL (Z, AIFCS, NAIF) - X*(0.25 +
     1  CSEVL (Z, AIGCS, NAIG)) )
      RETURN
C
 30   IF (X.GT.XMAX) GO TO 40
      AI = AIE(X) * EXP(-2.0*X*SQRT(X)/3.0)
      RETURN
C
 40   AI = 0.0
      CALL XERMSG ('SLATEC', 'AI', 'X SO BIG AI UNDERFLOWS', 1, 1)
      RETURN
C
      END
