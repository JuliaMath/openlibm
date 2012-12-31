*DECK R9ATN1
      FUNCTION R9ATN1 (X)
C***BEGIN PROLOGUE  R9ATN1
C***SUBSIDIARY
C***PURPOSE  Evaluate ATAN(X) from first order relative accuracy so that
C            ATAN(X) = X + X**3*R9ATN1(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      SINGLE PRECISION (R9ATN1-S, D9ATN1-D)
C***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB,
C             TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate  ATAN(X)  from first order, that is, evaluate
C (ATAN(X)-X)/X**3  with relative error accuracy so that
C        ATAN(X) = X + X**3*R9ATN1(X).
C
C Series for ATN1       on the interval  0.          to  1.00000D+00
C                                        with weighted error   2.21E-17
C                                         log weighted error  16.66
C                               significant figures required  15.44
C                                    decimal places required  17.32
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780401  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  R9ATN1
      DIMENSION ATN1CS(21)
      LOGICAL FIRST
      SAVE ATN1CS, NTATN1, XSML, XBIG, XMAX, FIRST
      DATA ATN1CS( 1) /   -.0328399753 5355202E0 /
      DATA ATN1CS( 2) /    .0583343234 3172412E0 /
      DATA ATN1CS( 3) /   -.0074003696 9671964E0 /
      DATA ATN1CS( 4) /    .0010097841 9933728E0 /
      DATA ATN1CS( 5) /   -.0001439787 1635652E0 /
      DATA ATN1CS( 6) /    .0000211451 2648992E0 /
      DATA ATN1CS( 7) /   -.0000031723 2107425E0 /
      DATA ATN1CS( 8) /    .0000004836 6203654E0 /
      DATA ATN1CS( 9) /   -.0000000746 7746546E0 /
      DATA ATN1CS(10) /    .0000000116 4800896E0 /
      DATA ATN1CS(11) /   -.0000000018 3208837E0 /
      DATA ATN1CS(12) /    .0000000002 9019082E0 /
      DATA ATN1CS(13) /   -.0000000000 4623885E0 /
      DATA ATN1CS(14) /    .0000000000 0740552E0 /
      DATA ATN1CS(15) /   -.0000000000 0119135E0 /
      DATA ATN1CS(16) /    .0000000000 0019240E0 /
      DATA ATN1CS(17) /   -.0000000000 0003118E0 /
      DATA ATN1CS(18) /    .0000000000 0000506E0 /
      DATA ATN1CS(19) /   -.0000000000 0000082E0 /
      DATA ATN1CS(20) /    .0000000000 0000013E0 /
      DATA ATN1CS(21) /   -.0000000000 0000002E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  R9ATN1
      IF (FIRST) THEN
         EPS = R1MACH(3)
         NTATN1 = INITS (ATN1CS, 21, 0.1*EPS)
C
         XSML = SQRT (0.1*EPS)
         XBIG = 1.571/SQRT(EPS)
         XMAX = 1.571/EPS
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.1.0) GO TO 20
C
      IF (Y.LE.XSML) R9ATN1 = -1.0/3.0
      IF (Y.LE.XSML) RETURN
C
      R9ATN1 = -0.25 + CSEVL (2.0*Y*Y-1., ATN1CS, NTATN1)
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'R9ATN1',
     +   'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG', 2, 2)
      IF (Y .GT. XBIG) CALL XERMSG ('SLATEC', 'R9ATN1',
     +   'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG', 1, 1)
C
      R9ATN1 = (ATAN(X) - X) / X**3
      RETURN
C
      END
