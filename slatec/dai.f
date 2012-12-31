*DECK DAI
      DOUBLE PRECISION FUNCTION DAI (X)
C***BEGIN PROLOGUE  DAI
C***PURPOSE  Evaluate the Airy function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10D
C***TYPE      DOUBLE PRECISION (AI-S, DAI-D)
C***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DAI(X) calculates the double precision Airy function for double
C precision argument X.
C
C Series for AIF        on the interval -1.00000E+00 to  1.00000E+00
C                                        with weighted error   8.37E-33
C                                         log weighted error  32.08
C                               significant figures required  30.87
C                                    decimal places required  32.63
C
C Series for AIG        on the interval -1.00000E+00 to  1.00000E+00
C                                        with weighted error   7.47E-34
C                                         log weighted error  33.13
C                               significant figures required  31.50
C                                    decimal places required  33.68
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9AIMP, DAIE, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  DAI
      DOUBLE PRECISION X, AIFCS(13), AIGCS(13), THETA, XM, XMAX, X3SML,
     1  Z, D1MACH, DCSEVL, DAIE, XMAXT
      LOGICAL FIRST
      SAVE AIFCS, AIGCS, NAIF, NAIG, X3SML, XMAX, FIRST
      DATA AIFCS(  1) / -.3797135849 6669997496 1970894694 14 D-1     /
      DATA AIFCS(  2) / +.5919188853 7263638574 3197280137 77 D-1     /
      DATA AIFCS(  3) / +.9862928057 7279975365 6038910440 60 D-3     /
      DATA AIFCS(  4) / +.6848843819 0765667554 8548301824 12 D-5     /
      DATA AIFCS(  5) / +.2594202596 2194713019 4892790814 03 D-7     /
      DATA AIFCS(  6) / +.6176612774 0813750329 4457496972 36 D-10    /
      DATA AIFCS(  7) / +.1009245417 2466117901 4295562246 01 D-12    /
      DATA AIFCS(  8) / +.1201479251 1179938141 2880332253 33 D-15    /
      DATA AIFCS(  9) / +.1088294558 8716991878 5252954666 66 D-18    /
      DATA AIFCS( 10) / +.7751377219 6684887039 2384000000 00 D-22    /
      DATA AIFCS( 11) / +.4454811203 7175638391 4666666666 66 D-25    /
      DATA AIFCS( 12) / +.2109284523 1692343466 6666666666 66 D-28    /
      DATA AIFCS( 13) / +.8370173591 0741333333 3333333333 33 D-32    /
      DATA AIGCS(  1) / +.1815236558 1161273011 5562099578 64 D-1     /
      DATA AIGCS(  2) / +.2157256316 6010755534 0306388199 68 D-1     /
      DATA AIGCS(  3) / +.2567835698 7483249659 0524280901 33 D-3     /
      DATA AIGCS(  4) / +.1426521411 9792403898 8294969217 21 D-5     /
      DATA AIGCS(  5) / +.4572114920 0180426070 4340975581 91 D-8     /
      DATA AIGCS(  6) / +.9525170843 5647098607 3922788405 92 D-11    /
      DATA AIGCS(  7) / +.1392563460 5771399051 1504206861 90 D-13    /
      DATA AIGCS(  8) / +.1507099914 2762379592 3069911386 66 D-16    /
      DATA AIGCS(  9) / +.1255914831 2567778822 7032053333 33 D-19    /
      DATA AIGCS( 10) / +.8306307377 0821340343 8293333333 33 D-23    /
      DATA AIGCS( 11) / +.4465753849 3718567445 3333333333 33 D-26    /
      DATA AIGCS( 12) / +.1990085503 4518869333 3333333333 33 D-29    /
      DATA AIGCS( 13) / +.7470288525 6533333333 3333333333 33 D-33    /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DAI
      IF (FIRST) THEN
         NAIF = INITDS (AIFCS, 13, 0.1*REAL(D1MACH(3)))
         NAIG = INITDS (AIGCS, 13, 0.1*REAL(D1MACH(3)))
C
         X3SML = D1MACH(3)**0.3334D0
         XMAXT = (-1.5D0*LOG(D1MACH(1)))**0.6667D0
         XMAX = XMAXT - XMAXT*LOG(XMAXT)/(4.0D0*SQRT(XMAXT)+1.0D0)
     *           - 0.01D0
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GE.(-1.D0)) GO TO 20
      CALL D9AIMP (X, XM, THETA)
      DAI = XM * COS(THETA)
      RETURN
C
 20   IF (X.GT.1.0D0) GO TO 30
      Z = 0.0D0
      IF (ABS(X).GT.X3SML) Z = X**3
      DAI = 0.375D0 + (DCSEVL (Z, AIFCS, NAIF) - X*(0.25D0 +
     1  DCSEVL (Z, AIGCS, NAIG)) )
      RETURN
C
 30   IF (X.GT.XMAX) GO TO 40
      DAI = DAIE(X) * EXP(-2.0D0*X*SQRT(X)/3.0D0)
      RETURN
C
 40   DAI = 0.0D0
      CALL XERMSG ('SLATEC', 'DAI', 'X SO BIG AI UNDERFLOWS', 1, 1)
      RETURN
C
      END
