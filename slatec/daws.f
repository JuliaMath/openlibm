*DECK DAWS
      FUNCTION DAWS (X)
C***BEGIN PROLOGUE  DAWS
C***PURPOSE  Compute Dawson's function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8C
C***TYPE      SINGLE PRECISION (DAWS-S, DDAWS-D)
C***KEYWORDS  DAWSON'S FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DAWS(X) calculates Dawson's integral for real argument X.
C
C Series for DAW        on the interval  0.          to  1.00000D+00
C                                        with weighted error   3.83E-17
C                                         log weighted error  16.42
C                               significant figures required  15.78
C                                    decimal places required  16.97
C
C Series for DAW2       on the interval  0.          to  1.60000D+01
C                                        with weighted error   5.17E-17
C                                         log weighted error  16.29
C                               significant figures required  15.90
C                                    decimal places required  17.02
C
C Series for DAWA       on the interval  0.          to  6.25000D-02
C                                        with weighted error   2.24E-17
C                                         log weighted error  16.65
C                               significant figures required  14.73
C                                    decimal places required  17.36
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  DAWS
      DIMENSION DAWCS(13), DAW2CS(29), DAWACS(26)
      LOGICAL FIRST
      SAVE DAWCS, DAW2CS, DAWACS, NTDAW, NTDAW2, NTDAWA,
     1 XSML, XBIG, XMAX, FIRST
      DATA DAWCS( 1) /   -.0063517343 75145949E0 /
      DATA DAWCS( 2) /   -.2294071479 6773869E0 /
      DATA DAWCS( 3) /    .0221305009 39084764E0 /
      DATA DAWCS( 4) /   -.0015492654 53892985E0 /
      DATA DAWCS( 5) /    .0000849732 77156849E0 /
      DATA DAWCS( 6) /   -.0000038282 66270972E0 /
      DATA DAWCS( 7) /    .0000001462 85480625E0 /
      DATA DAWCS( 8) /   -.0000000048 51982381E0 /
      DATA DAWCS( 9) /    .0000000001 42146357E0 /
      DATA DAWCS(10) /   -.0000000000 03728836E0 /
      DATA DAWCS(11) /    .0000000000 00088549E0 /
      DATA DAWCS(12) /   -.0000000000 00001920E0 /
      DATA DAWCS(13) /    .0000000000 00000038E0 /
      DATA DAW2CS( 1) /   -.0568865441 05215527E0 /
      DATA DAW2CS( 2) /   -.3181134699 6168131E0 /
      DATA DAW2CS( 3) /    .2087384541 3642237E0 /
      DATA DAW2CS( 4) /   -.1247540991 3779131E0 /
      DATA DAW2CS( 5) /    .0678693051 86676777E0 /
      DATA DAW2CS( 6) /   -.0336591448 95270940E0 /
      DATA DAW2CS( 7) /    .0152607812 71987972E0 /
      DATA DAW2CS( 8) /   -.0063483709 62596214E0 /
      DATA DAW2CS( 9) /    .0024326740 92074852E0 /
      DATA DAW2CS(10) /   -.0008621954 14910650E0 /
      DATA DAW2CS(11) /    .0002837657 33363216E0 /
      DATA DAW2CS(12) /   -.0000870575 49874170E0 /
      DATA DAW2CS(13) /    .0000249868 49985481E0 /
      DATA DAW2CS(14) /   -.0000067319 28676416E0 /
      DATA DAW2CS(15) /    .0000017078 57878557E0 /
      DATA DAW2CS(16) /   -.0000004091 75512264E0 /
      DATA DAW2CS(17) /    .0000000928 28292216E0 /
      DATA DAW2CS(18) /   -.0000000199 91403610E0 /
      DATA DAW2CS(19) /    .0000000040 96349064E0 /
      DATA DAW2CS(20) /   -.0000000008 00324095E0 /
      DATA DAW2CS(21) /    .0000000001 49385031E0 /
      DATA DAW2CS(22) /   -.0000000000 26687999E0 /
      DATA DAW2CS(23) /    .0000000000 04571221E0 /
      DATA DAW2CS(24) /   -.0000000000 00751873E0 /
      DATA DAW2CS(25) /    .0000000000 00118931E0 /
      DATA DAW2CS(26) /   -.0000000000 00018116E0 /
      DATA DAW2CS(27) /    .0000000000 00002661E0 /
      DATA DAW2CS(28) /   -.0000000000 00000377E0 /
      DATA DAW2CS(29) /    .0000000000 00000051E0 /
      DATA DAWACS( 1) /    .0169048563 7765704E0 /
      DATA DAWACS( 2) /    .0086832522 7840695E0 /
      DATA DAWACS( 3) /    .0002424864 0424177E0 /
      DATA DAWACS( 4) /    .0000126118 2399572E0 /
      DATA DAWACS( 5) /    .0000010664 5331463E0 /
      DATA DAWACS( 6) /    .0000001358 1597947E0 /
      DATA DAWACS( 7) /    .0000000217 1042356E0 /
      DATA DAWACS( 8) /    .0000000028 6701050E0 /
      DATA DAWACS( 9) /   -.0000000001 9013363E0 /
      DATA DAWACS(10) /   -.0000000003 0977804E0 /
      DATA DAWACS(11) /   -.0000000001 0294148E0 /
      DATA DAWACS(12) /   -.0000000000 0626035E0 /
      DATA DAWACS(13) /    .0000000000 0856313E0 /
      DATA DAWACS(14) /    .0000000000 0303304E0 /
      DATA DAWACS(15) /   -.0000000000 0025236E0 /
      DATA DAWACS(16) /   -.0000000000 0042106E0 /
      DATA DAWACS(17) /   -.0000000000 0004431E0 /
      DATA DAWACS(18) /    .0000000000 0004911E0 /
      DATA DAWACS(19) /    .0000000000 0001235E0 /
      DATA DAWACS(20) /   -.0000000000 0000578E0 /
      DATA DAWACS(21) /   -.0000000000 0000228E0 /
      DATA DAWACS(22) /    .0000000000 0000076E0 /
      DATA DAWACS(23) /    .0000000000 0000038E0 /
      DATA DAWACS(24) /   -.0000000000 0000011E0 /
      DATA DAWACS(25) /   -.0000000000 0000006E0 /
      DATA DAWACS(26) /    .0000000000 0000002E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DAWS
      IF (FIRST) THEN
         EPS = R1MACH(3)
         NTDAW  = INITS (DAWCS,  13, 0.1*EPS)
         NTDAW2 = INITS (DAW2CS, 29, 0.1*EPS)
         NTDAWA = INITS (DAWACS, 26, 0.1*EPS)
C
         XSML = SQRT (1.5*EPS)
         XBIG = SQRT (0.5/EPS)
         XMAX = EXP (MIN (-LOG(2.*R1MACH(1)), LOG(R1MACH(2))) - 1.0)
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.1.0) GO TO 20
C
      DAWS = X
      IF (Y.LE.XSML) RETURN
C
      DAWS = X * (0.75 + CSEVL (2.0*Y*Y-1.0, DAWCS, NTDAW))
      RETURN
C
 20   IF (Y.GT.4.0) GO TO 30
      DAWS = X * (0.25 + CSEVL (0.125*Y*Y-1.0, DAW2CS, NTDAW2))
      RETURN
C
 30   IF (Y.GT.XMAX) GO TO 40
      DAWS = 0.5/X
      IF (Y.GT.XBIG) RETURN
C
      DAWS = (0.5 + CSEVL (32.0/Y**2-1.0, DAWACS, NTDAWA)) / X
      RETURN
C
 40   CALL XERMSG ('SLATEC', 'DAWS', 'ABS(X) SO LARGE DAWS UNDERFLOWS',
     +   1, 1)
      DAWS = 0.0
      RETURN
C
      END
