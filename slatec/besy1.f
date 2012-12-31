*DECK BESY1
      FUNCTION BESY1 (X)
C***BEGIN PROLOGUE  BESY1
C***PURPOSE  Compute the Bessel function of the second kind of order
C            one.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      SINGLE PRECISION (BESY1-S, DBESY1-D)
C***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ONE, SECOND KIND,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESY1(X) calculates the Bessel function of the second kind of
C order one for real argument X.
C
C Series for BY1        on the interval  0.          to  1.60000D+01
C                                        with weighted error   1.87E-18
C                                         log weighted error  17.73
C                               significant figures required  17.83
C                                    decimal places required  18.30
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
C***ROUTINES CALLED  BESJ1, CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C***END PROLOGUE  BESY1
      DIMENSION BY1CS(14), BM1CS(21), BTH1CS(24)
      LOGICAL FIRST
      SAVE BY1CS, BM1CS, BTH1CS, TWODPI, PI4,
     1 NTY1, NTM1, NTTH1, XMIN, XSML, XMAX, FIRST
      DATA BY1CS( 1) /    .0320804710 0611908629E0 /
      DATA BY1CS( 2) /   1.2627078974 33500450E0 /
      DATA BY1CS( 3) /    .0064999618 9992317500E0 /
      DATA BY1CS( 4) /   -.0893616452 8860504117E0 /
      DATA BY1CS( 5) /    .0132508812 2175709545E0 /
      DATA BY1CS( 6) /   -.0008979059 1196483523E0 /
      DATA BY1CS( 7) /    .0000364736 1487958306E0 /
      DATA BY1CS( 8) /   -.0000010013 7438166600E0 /
      DATA BY1CS( 9) /    .0000000199 4539657390E0 /
      DATA BY1CS(10) /   -.0000000003 0230656018E0 /
      DATA BY1CS(11) /    .0000000000 0360987815E0 /
      DATA BY1CS(12) /   -.0000000000 0003487488E0 /
      DATA BY1CS(13) /    .0000000000 0000027838E0 /
      DATA BY1CS(14) /   -.0000000000 0000000186E0 /
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
      DATA TWODPI / 0.6366197723 6758134E0 /
      DATA PI4 / 0.7853981633 9744831E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESY1
      IF (FIRST) THEN
         NTY1 = INITS (BY1CS, 14, 0.1*R1MACH(3))
         NTM1 = INITS (BM1CS, 21, 0.1*R1MACH(3))
         NTTH1 = INITS (BTH1CS, 24, 0.1*R1MACH(3))
C
         XMIN = 1.571*EXP ( MAX(LOG(R1MACH(1)), -LOG(R1MACH(2)))+.01)
         XSML = SQRT (4.0*R1MACH(3))
         XMAX = 1.0/R1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. 0.) CALL XERMSG ('SLATEC', 'BESY1',
     +   'X IS ZERO OR NEGATIVE', 1, 2)
      IF (X.GT.4.0) GO TO 20
C
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'BESY1',
     +   'X SO SMALL Y1 OVERFLOWS', 3, 2)
      Y = 0.
      IF (X.GT.XSML) Y = X*X
      BESY1 = TWODPI*LOG(0.5*X)*BESJ1(X) +
     1  (0.5 + CSEVL (.125*Y-1., BY1CS, NTY1))/X
      RETURN
C
 20   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'BESY1',
     +   'NO PRECISION BECAUSE X IS BIG', 2, 2)
C
      Z = 32.0/X**2 - 1.0
      AMPL = (0.75 + CSEVL (Z, BM1CS, NTM1)) / SQRT(X)
      THETA = X - 3.0*PI4 + CSEVL (Z, BTH1CS, NTTH1) / X
      BESY1 = AMPL * SIN (THETA)
C
      RETURN
      END
