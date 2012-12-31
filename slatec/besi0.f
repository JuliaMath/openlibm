*DECK BESI0
      FUNCTION BESI0 (X)
C***BEGIN PROLOGUE  BESI0
C***PURPOSE  Compute the hyperbolic Bessel function of the first kind
C            of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESI0-S, DBESI0-D)
C***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
C             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESI0(X) computes the modified (hyperbolic) Bessel function
C of the first kind of order zero and real argument X.
C
C Series for BI0        on the interval  0.          to  9.00000D+00
C                                        with weighted error   2.46E-18
C                                         log weighted error  17.61
C                               significant figures required  17.90
C                                    decimal places required  18.15
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  BESI0E, CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESI0
      DIMENSION BI0CS(12)
      LOGICAL FIRST
      SAVE BI0CS, NTI0, XSML, XMAX, FIRST
      DATA BI0CS( 1) /   -.0766054725 2839144951E0 /
      DATA BI0CS( 2) /   1.9273379539 93808270E0 /
      DATA BI0CS( 3) /    .2282644586 920301339E0 /
      DATA BI0CS( 4) /    .0130489146 6707290428E0 /
      DATA BI0CS( 5) /    .0004344270 9008164874E0 /
      DATA BI0CS( 6) /    .0000094226 5768600193E0 /
      DATA BI0CS( 7) /    .0000001434 0062895106E0 /
      DATA BI0CS( 8) /    .0000000016 1384906966E0 /
      DATA BI0CS( 9) /    .0000000000 1396650044E0 /
      DATA BI0CS(10) /    .0000000000 0009579451E0 /
      DATA BI0CS(11) /    .0000000000 0000053339E0 /
      DATA BI0CS(12) /    .0000000000 0000000245E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESI0
      IF (FIRST) THEN
         NTI0 = INITS (BI0CS, 12, 0.1*R1MACH(3))
         XSML = SQRT (4.5*R1MACH(3))
         XMAX = LOG (R1MACH(2))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0) GO TO 20
C
      BESI0 = 1.0
      IF (Y.GT.XSML) BESI0 = 2.75 + CSEVL (Y*Y/4.5-1.0, BI0CS, NTI0)
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'BESI0',
     +   'ABS(X) SO BIG I0 OVERFLOWS', 1, 2)
C
      BESI0 = EXP(Y) * BESI0E(X)
C
      RETURN
      END
