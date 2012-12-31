*DECK BESJ0
      FUNCTION BESJ0 (X)
C***BEGIN PROLOGUE  BESJ0
C***PURPOSE  Compute the Bessel function of the first kind of order
C            zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      SINGLE PRECISION (BESJ0-S, DBESJ0-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESJ0(X) calculates the Bessel function of the first kind of
C order zero for real argument X.
C
C Series for BJ0        on the interval  0.          to  1.60000D+01
C                                        with weighted error   7.47E-18
C                                         log weighted error  17.13
C                               significant figures required  16.98
C                                    decimal places required  17.68
C
C Series for BM0        on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.98E-17
C                                         log weighted error  16.30
C                               significant figures required  14.97
C                                    decimal places required  16.96
C
C Series for BTH0       on the interval  0.          to  6.25000D-02
C                                        with weighted error   3.67E-17
C                                         log weighted error  16.44
C                               significant figures required  15.53
C                                    decimal places required  17.13
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890210  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESJ0
      DIMENSION BJ0CS(13), BM0CS(21), BTH0CS(24)
      LOGICAL FIRST
      SAVE BJ0CS, BM0CS, BTH0CS, PI4, NTJ0, NTM0, NTTH0, XSML, XMAX,
     1   FIRST
      DATA BJ0CS( 1) /    .1002541619 68939137E0 /
      DATA BJ0CS( 2) /   -.6652230077 64405132E0 /
      DATA BJ0CS( 3) /    .2489837034 98281314E0 /
      DATA BJ0CS( 4) /   -.0332527231 700357697E0 /
      DATA BJ0CS( 5) /    .0023114179 304694015E0 /
      DATA BJ0CS( 6) /   -.0000991127 741995080E0 /
      DATA BJ0CS( 7) /    .0000028916 708643998E0 /
      DATA BJ0CS( 8) /   -.0000000612 108586630E0 /
      DATA BJ0CS( 9) /    .0000000009 838650793E0 /
      DATA BJ0CS(10) /   -.0000000000 124235515E0 /
      DATA BJ0CS(11) /    .0000000000 001265433E0 /
      DATA BJ0CS(12) /   -.0000000000 000010619E0 /
      DATA BJ0CS(13) /    .0000000000 000000074E0 /
      DATA BM0CS( 1) /    .0928496163 7381644E0 /
      DATA BM0CS( 2) /   -.0014298770 7403484E0 /
      DATA BM0CS( 3) /    .0000283057 9271257E0 /
      DATA BM0CS( 4) /   -.0000014330 0611424E0 /
      DATA BM0CS( 5) /    .0000001202 8628046E0 /
      DATA BM0CS( 6) /   -.0000000139 7113013E0 /
      DATA BM0CS( 7) /    .0000000020 4076188E0 /
      DATA BM0CS( 8) /   -.0000000003 5399669E0 /
      DATA BM0CS( 9) /    .0000000000 7024759E0 /
      DATA BM0CS(10) /   -.0000000000 1554107E0 /
      DATA BM0CS(11) /    .0000000000 0376226E0 /
      DATA BM0CS(12) /   -.0000000000 0098282E0 /
      DATA BM0CS(13) /    .0000000000 0027408E0 /
      DATA BM0CS(14) /   -.0000000000 0008091E0 /
      DATA BM0CS(15) /    .0000000000 0002511E0 /
      DATA BM0CS(16) /   -.0000000000 0000814E0 /
      DATA BM0CS(17) /    .0000000000 0000275E0 /
      DATA BM0CS(18) /   -.0000000000 0000096E0 /
      DATA BM0CS(19) /    .0000000000 0000034E0 /
      DATA BM0CS(20) /   -.0000000000 0000012E0 /
      DATA BM0CS(21) /    .0000000000 0000004E0 /
      DATA BTH0CS( 1) /   -.2463916377 4300119E0 /
      DATA BTH0CS( 2) /    .0017370983 07508963E0 /
      DATA BTH0CS( 3) /   -.0000621836 33402968E0 /
      DATA BTH0CS( 4) /    .0000043680 50165742E0 /
      DATA BTH0CS( 5) /   -.0000004560 93019869E0 /
      DATA BTH0CS( 6) /    .0000000621 97400101E0 /
      DATA BTH0CS( 7) /   -.0000000103 00442889E0 /
      DATA BTH0CS( 8) /    .0000000019 79526776E0 /
      DATA BTH0CS( 9) /   -.0000000004 28198396E0 /
      DATA BTH0CS(10) /    .0000000001 02035840E0 /
      DATA BTH0CS(11) /   -.0000000000 26363898E0 /
      DATA BTH0CS(12) /    .0000000000 07297935E0 /
      DATA BTH0CS(13) /   -.0000000000 02144188E0 /
      DATA BTH0CS(14) /    .0000000000 00663693E0 /
      DATA BTH0CS(15) /   -.0000000000 00215126E0 /
      DATA BTH0CS(16) /    .0000000000 00072659E0 /
      DATA BTH0CS(17) /   -.0000000000 00025465E0 /
      DATA BTH0CS(18) /    .0000000000 00009229E0 /
      DATA BTH0CS(19) /   -.0000000000 00003448E0 /
      DATA BTH0CS(20) /    .0000000000 00001325E0 /
      DATA BTH0CS(21) /   -.0000000000 00000522E0 /
      DATA BTH0CS(22) /    .0000000000 00000210E0 /
      DATA BTH0CS(23) /   -.0000000000 00000087E0 /
      DATA BTH0CS(24) /    .0000000000 00000036E0 /
      DATA PI4 / 0.7853981633 9744831E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESJ0
      IF (FIRST) THEN
         NTJ0 = INITS (BJ0CS, 13, 0.1*R1MACH(3))
         NTM0 = INITS (BM0CS, 21, 0.1*R1MACH(3))
         NTTH0 = INITS (BTH0CS, 24, 0.1*R1MACH(3))
C
         XSML = SQRT (8.0*R1MACH(3))
         XMAX = 1.0/R1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.4.0) GO TO 20
C
      BESJ0 = 1.0
      IF (Y.GT.XSML) BESJ0 = CSEVL (.125*Y*Y-1., BJ0CS, NTJ0)
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'BESJ0',
     +   'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 1, 2)
C
      Z = 32.0/Y**2 - 1.0
      AMPL = (0.75 + CSEVL (Z, BM0CS, NTM0)) / SQRT(Y)
      THETA = Y - PI4 + CSEVL (Z, BTH0CS, NTTH0) / Y
      BESJ0 = AMPL * COS (THETA)
C
      RETURN
      END
