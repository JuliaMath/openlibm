*DECK BESJ1
      FUNCTION BESJ1 (X)
C***BEGIN PROLOGUE  BESJ1
C***PURPOSE  Compute the Bessel function of the first kind of order one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      SINGLE PRECISION (BESJ1-S, DBESJ1-D)
C***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESJ1(X) calculates the Bessel function of the first kind of
C order one for real argument X.
C
C Series for BJ1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   4.48E-17
C                                         log weighted error  16.35
C                               significant figures required  15.77
C                                    decimal places required  16.89
C
C Series for BM1        on the interval  0.          to  6.25000D-02
C                                        with weighted error   5.61E-17
C                                         log weighted error  16.25
C                               significant figures required  14.97
C                                    decimal places required  16.91
C
C Series for BTH1       on the interval  0.          to  6.25000D-02
C                                        with weighted error   4.10E-17
C                                         log weighted error  16.39
C                               significant figures required  15.96
C                                    decimal places required  17.08
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   780601  DATE WRITTEN
C   890210  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESJ1
      DIMENSION BJ1CS(12), BM1CS(21), BTH1CS(24)
      LOGICAL FIRST
      SAVE BJ1CS, BM1CS, BTH1CS, PI4, NTJ1, NTM1, NTTH1,
     1 XSML, XMIN, XMAX, FIRST
      DATA BJ1CS( 1) /   -.1172614151 3332787E0 /
      DATA BJ1CS( 2) /   -.2536152183 0790640E0 /
      DATA BJ1CS( 3) /    .0501270809 84469569E0 /
      DATA BJ1CS( 4) /   -.0046315148 09625081E0 /
      DATA BJ1CS( 5) /    .0002479962 29415914E0 /
      DATA BJ1CS( 6) /   -.0000086789 48686278E0 /
      DATA BJ1CS( 7) /    .0000002142 93917143E0 /
      DATA BJ1CS( 8) /   -.0000000039 36093079E0 /
      DATA BJ1CS( 9) /    .0000000000 55911823E0 /
      DATA BJ1CS(10) /   -.0000000000 00632761E0 /
      DATA BJ1CS(11) /    .0000000000 00005840E0 /
      DATA BJ1CS(12) /   -.0000000000 00000044E0 /
      DATA BM1CS( 1) /    .1047362510 931285E0 /
      DATA BM1CS( 2) /    .0044244389 3702345E0 /
      DATA BM1CS( 3) /   -.0000566163 9504035E0 /
      DATA BM1CS( 4) /    .0000023134 9417339E0 /
      DATA BM1CS( 5) /   -.0000001737 7182007E0 /
      DATA BM1CS( 6) /    .0000000189 3209930E0 /
      DATA BM1CS( 7) /   -.0000000026 5416023E0 /
      DATA BM1CS( 8) /    .0000000004 4740209E0 /
      DATA BM1CS( 9) /   -.0000000000 8691795E0 /
      DATA BM1CS(10) /    .0000000000 1891492E0 /
      DATA BM1CS(11) /   -.0000000000 0451884E0 /
      DATA BM1CS(12) /    .0000000000 0116765E0 /
      DATA BM1CS(13) /   -.0000000000 0032265E0 /
      DATA BM1CS(14) /    .0000000000 0009450E0 /
      DATA BM1CS(15) /   -.0000000000 0002913E0 /
      DATA BM1CS(16) /    .0000000000 0000939E0 /
      DATA BM1CS(17) /   -.0000000000 0000315E0 /
      DATA BM1CS(18) /    .0000000000 0000109E0 /
      DATA BM1CS(19) /   -.0000000000 0000039E0 /
      DATA BM1CS(20) /    .0000000000 0000014E0 /
      DATA BM1CS(21) /   -.0000000000 0000005E0 /
      DATA BTH1CS( 1) /    .7406014102 6313850E0 /
      DATA BTH1CS( 2) /   -.0045717556 59637690E0 /
      DATA BTH1CS( 3) /    .0001198185 10964326E0 /
      DATA BTH1CS( 4) /   -.0000069645 61891648E0 /
      DATA BTH1CS( 5) /    .0000006554 95621447E0 /
      DATA BTH1CS( 6) /   -.0000000840 66228945E0 /
      DATA BTH1CS( 7) /    .0000000133 76886564E0 /
      DATA BTH1CS( 8) /   -.0000000024 99565654E0 /
      DATA BTH1CS( 9) /    .0000000005 29495100E0 /
      DATA BTH1CS(10) /   -.0000000001 24135944E0 /
      DATA BTH1CS(11) /    .0000000000 31656485E0 /
      DATA BTH1CS(12) /   -.0000000000 08668640E0 /
      DATA BTH1CS(13) /    .0000000000 02523758E0 /
      DATA BTH1CS(14) /   -.0000000000 00775085E0 /
      DATA BTH1CS(15) /    .0000000000 00249527E0 /
      DATA BTH1CS(16) /   -.0000000000 00083773E0 /
      DATA BTH1CS(17) /    .0000000000 00029205E0 /
      DATA BTH1CS(18) /   -.0000000000 00010534E0 /
      DATA BTH1CS(19) /    .0000000000 00003919E0 /
      DATA BTH1CS(20) /   -.0000000000 00001500E0 /
      DATA BTH1CS(21) /    .0000000000 00000589E0 /
      DATA BTH1CS(22) /   -.0000000000 00000237E0 /
      DATA BTH1CS(23) /    .0000000000 00000097E0 /
      DATA BTH1CS(24) /   -.0000000000 00000040E0 /
      DATA PI4 / 0.7853981633 9744831E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESJ1
      IF (FIRST) THEN
         NTJ1 = INITS (BJ1CS, 12, 0.1*R1MACH(3))
         NTM1 = INITS (BM1CS, 21, 0.1*R1MACH(3))
         NTTH1 = INITS (BTH1CS, 24, 0.1*R1MACH(3))
C
         XSML = SQRT (8.0*R1MACH(3))
         XMIN = 2.0*R1MACH(1)
         XMAX = 1.0/R1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.4.0) GO TO 20
C
      BESJ1 = 0.
      IF (Y.EQ.0.0) RETURN
      IF (Y .LE. XMIN) CALL XERMSG ('SLATEC', 'BESJ1',
     +   'ABS(X) SO SMALL J1 UNDERFLOWS', 1, 1)
      IF (Y.GT.XMIN) BESJ1 = 0.5*X
      IF (Y.GT.XSML) BESJ1 = X * (.25 + CSEVL(.125*Y*Y-1., BJ1CS, NTJ1))
      RETURN
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'BESJ1',
     +   'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 2, 2)
      Z = 32.0/Y**2 - 1.0
      AMPL = (0.75 + CSEVL (Z, BM1CS, NTM1)) / SQRT(Y)
      THETA = Y - 3.0*PI4 + CSEVL (Z, BTH1CS, NTTH1) / Y
      BESJ1 = SIGN (AMPL, X) * COS (THETA)
C
      RETURN
      END
