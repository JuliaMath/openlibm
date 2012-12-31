*DECK E1
      FUNCTION E1 (X)
C***BEGIN PROLOGUE  E1
C***PURPOSE  Compute the exponential integral E1(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C5
C***TYPE      SINGLE PRECISION (E1-S, DE1-D)
C***KEYWORDS  E1 FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C E1 calculates the single precision exponential integral, E1(X), for
C positive single precision argument X and the Cauchy principal value
C for negative X.  If principal values are used everywhere, then, for
C all X,
C
C    E1(X) = -Ei(-X)
C or
C    Ei(X) = -E1(-X).
C
C
C Series for AE11       on the interval -1.00000D-01 to  0.
C                                        with weighted error   1.76E-17
C                                         log weighted error  16.75
C                               significant figures required  15.70
C                                    decimal places required  17.55
C
C
C Series for AE12       on the interval -2.50000D-01 to -1.00000D-01
C                                        with weighted error   5.83E-17
C                                         log weighted error  16.23
C                               significant figures required  15.76
C                                    decimal places required  16.93
C
C
C Series for E11        on the interval -4.00000D+00 to -1.00000D+00
C                                        with weighted error   1.08E-18
C                                         log weighted error  17.97
C                               significant figures required  19.02
C                                    decimal places required  18.61
C
C
C Series for E12        on the interval -1.00000D+00 to  1.00000D+00
C                                        with weighted error   3.15E-18
C                                         log weighted error  17.50
C                        approx significant figures required  15.8
C                                    decimal places required  18.10
C
C
C Series for AE13       on the interval  2.50000D-01 to  1.00000D+00
C                                        with weighted error   2.34E-17
C                                         log weighted error  16.63
C                               significant figures required  16.14
C                                    decimal places required  17.33
C
C
C Series for AE14       on the interval  0.          to  2.50000D-01
C                                        with weighted error   5.41E-17
C                                         log weighted error  16.27
C                               significant figures required  15.38
C                                    decimal places required  16.97
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891115  Modified prologue description.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  E1
      DIMENSION AE11CS(39), AE12CS(25), E11CS(19), E12CS(16),
     1  AE13CS(25), AE14CS(26)
      LOGICAL FIRST
      SAVE AE11CS, AE12CS, E11CS, E12CS, AE13CS, AE14CS,
     1 NTAE11, NTAE12, NTE11, NTE12, NTAE13, NTAE14, XMAX, FIRST
      DATA AE11CS( 1) /    .1215032397 1606579E0 /
      DATA AE11CS( 2) /   -.0650887785 13550150E0 /
      DATA AE11CS( 3) /    .0048976513 57459670E0 /
      DATA AE11CS( 4) /   -.0006492378 43027216E0 /
      DATA AE11CS( 5) /    .0000938404 34587471E0 /
      DATA AE11CS( 6) /    .0000004202 36380882E0 /
      DATA AE11CS( 7) /   -.0000081133 74735904E0 /
      DATA AE11CS( 8) /    .0000028042 47688663E0 /
      DATA AE11CS( 9) /    .0000000564 87164441E0 /
      DATA AE11CS(10) /   -.0000003448 09174450E0 /
      DATA AE11CS(11) /    .0000000582 09273578E0 /
      DATA AE11CS(12) /    .0000000387 11426349E0 /
      DATA AE11CS(13) /   -.0000000124 53235014E0 /
      DATA AE11CS(14) /   -.0000000051 18504888E0 /
      DATA AE11CS(15) /    .0000000021 48771527E0 /
      DATA AE11CS(16) /    .0000000008 68459898E0 /
      DATA AE11CS(17) /   -.0000000003 43650105E0 /
      DATA AE11CS(18) /   -.0000000001 79796603E0 /
      DATA AE11CS(19) /    .0000000000 47442060E0 /
      DATA AE11CS(20) /    .0000000000 40423282E0 /
      DATA AE11CS(21) /   -.0000000000 03543928E0 /
      DATA AE11CS(22) /   -.0000000000 08853444E0 /
      DATA AE11CS(23) /   -.0000000000 00960151E0 /
      DATA AE11CS(24) /    .0000000000 01692921E0 /
      DATA AE11CS(25) /    .0000000000 00607990E0 /
      DATA AE11CS(26) /   -.0000000000 00224338E0 /
      DATA AE11CS(27) /   -.0000000000 00200327E0 /
      DATA AE11CS(28) /   -.0000000000 00006246E0 /
      DATA AE11CS(29) /    .0000000000 00045571E0 /
      DATA AE11CS(30) /    .0000000000 00016383E0 /
      DATA AE11CS(31) /   -.0000000000 00005561E0 /
      DATA AE11CS(32) /   -.0000000000 00006074E0 /
      DATA AE11CS(33) /   -.0000000000 00000862E0 /
      DATA AE11CS(34) /    .0000000000 00001223E0 /
      DATA AE11CS(35) /    .0000000000 00000716E0 /
      DATA AE11CS(36) /   -.0000000000 00000024E0 /
      DATA AE11CS(37) /   -.0000000000 00000201E0 /
      DATA AE11CS(38) /   -.0000000000 00000082E0 /
      DATA AE11CS(39) /    .0000000000 00000017E0 /
      DATA AE12CS( 1) /    .5824174951 3472674E0 /
      DATA AE12CS( 2) /   -.1583488509 0578275E0 /
      DATA AE12CS( 3) /   -.0067642755 90323141E0 /
      DATA AE12CS( 4) /    .0051258439 50185725E0 /
      DATA AE12CS( 5) /    .0004352324 92169391E0 /
      DATA AE12CS( 6) /   -.0001436133 66305483E0 /
      DATA AE12CS( 7) /   -.0000418013 20556301E0 /
      DATA AE12CS( 8) /   -.0000027133 95758640E0 /
      DATA AE12CS( 9) /    .0000011513 81913647E0 /
      DATA AE12CS(10) /    .0000004206 50022012E0 /
      DATA AE12CS(11) /    .0000000665 81901391E0 /
      DATA AE12CS(12) /    .0000000006 62143777E0 /
      DATA AE12CS(13) /   -.0000000028 44104870E0 /
      DATA AE12CS(14) /   -.0000000009 40724197E0 /
      DATA AE12CS(15) /   -.0000000001 77476602E0 /
      DATA AE12CS(16) /   -.0000000000 15830222E0 /
      DATA AE12CS(17) /    .0000000000 02905732E0 /
      DATA AE12CS(18) /    .0000000000 01769356E0 /
      DATA AE12CS(19) /    .0000000000 00492735E0 /
      DATA AE12CS(20) /    .0000000000 00093709E0 /
      DATA AE12CS(21) /    .0000000000 00010707E0 /
      DATA AE12CS(22) /   -.0000000000 00000537E0 /
      DATA AE12CS(23) /   -.0000000000 00000716E0 /
      DATA AE12CS(24) /   -.0000000000 00000244E0 /
      DATA AE12CS(25) /   -.0000000000 00000058E0 /
      DATA E11CS( 1) / -16.1134616555 71494026E0 /
      DATA E11CS( 2) /   7.7940727787 426802769E0 /
      DATA E11CS( 3) /  -1.9554058188 631419507E0 /
      DATA E11CS( 4) /    .3733729386 6277945612E0 /
      DATA E11CS( 5) /   -.0569250319 1092901938E0 /
      DATA E11CS( 6) /    .0072110777 6966009185E0 /
      DATA E11CS( 7) /   -.0007810490 1449841593E0 /
      DATA E11CS( 8) /    .0000738809 3356262168E0 /
      DATA E11CS( 9) /   -.0000062028 6187580820E0 /
      DATA E11CS(10) /    .0000004681 6002303176E0 /
      DATA E11CS(11) /   -.0000000320 9288853329E0 /
      DATA E11CS(12) /    .0000000020 1519974874E0 /
      DATA E11CS(13) /   -.0000000001 1673686816E0 /
      DATA E11CS(14) /    .0000000000 0627627066E0 /
      DATA E11CS(15) /   -.0000000000 0031481541E0 /
      DATA E11CS(16) /    .0000000000 0001479904E0 /
      DATA E11CS(17) /   -.0000000000 0000065457E0 /
      DATA E11CS(18) /    .0000000000 0000002733E0 /
      DATA E11CS(19) /   -.0000000000 0000000108E0 /
      DATA E12CS( 1) /  -0.0373902147 92202795E0 /
      DATA E12CS( 2) /   0.0427239860 62209577E0 /
      DATA E12CS( 3) /   -.1303182079 849700544E0 /
      DATA E12CS( 4) /    .0144191240 2469889073E0 /
      DATA E12CS( 5) /   -.0013461707 8051068022E0 /
      DATA E12CS( 6) /    .0001073102 9253063780E0 /
      DATA E12CS( 7) /   -.0000074299 9951611943E0 /
      DATA E12CS( 8) /    .0000004537 7325690753E0 /
      DATA E12CS( 9) /   -.0000000247 6417211390E0 /
      DATA E12CS(10) /    .0000000012 2076581374E0 /
      DATA E12CS(11) /   -.0000000000 5485141480E0 /
      DATA E12CS(12) /    .0000000000 0226362142E0 /
      DATA E12CS(13) /   -.0000000000 0008635897E0 /
      DATA E12CS(14) /    .0000000000 0000306291E0 /
      DATA E12CS(15) /   -.0000000000 0000010148E0 /
      DATA E12CS(16) /    .0000000000 0000000315E0 /
      DATA AE13CS( 1) /   -.6057732466 4060346E0 /
      DATA AE13CS( 2) /   -.1125352434 8366090E0 /
      DATA AE13CS( 3) /    .0134322662 47902779E0 /
      DATA AE13CS( 4) /   -.0019268451 87381145E0 /
      DATA AE13CS( 5) /    .0003091183 37720603E0 /
      DATA AE13CS( 6) /   -.0000535641 32129618E0 /
      DATA AE13CS( 7) /    .0000098278 12880247E0 /
      DATA AE13CS( 8) /   -.0000018853 68984916E0 /
      DATA AE13CS( 9) /    .0000003749 43193568E0 /
      DATA AE13CS(10) /   -.0000000768 23455870E0 /
      DATA AE13CS(11) /    .0000000161 43270567E0 /
      DATA AE13CS(12) /   -.0000000034 66802211E0 /
      DATA AE13CS(13) /    .0000000007 58754209E0 /
      DATA AE13CS(14) /   -.0000000001 68864333E0 /
      DATA AE13CS(15) /    .0000000000 38145706E0 /
      DATA AE13CS(16) /   -.0000000000 08733026E0 /
      DATA AE13CS(17) /    .0000000000 02023672E0 /
      DATA AE13CS(18) /   -.0000000000 00474132E0 /
      DATA AE13CS(19) /    .0000000000 00112211E0 /
      DATA AE13CS(20) /   -.0000000000 00026804E0 /
      DATA AE13CS(21) /    .0000000000 00006457E0 /
      DATA AE13CS(22) /   -.0000000000 00001568E0 /
      DATA AE13CS(23) /    .0000000000 00000383E0 /
      DATA AE13CS(24) /   -.0000000000 00000094E0 /
      DATA AE13CS(25) /    .0000000000 00000023E0 /
      DATA AE14CS( 1) /   -.1892918000 753017E0 /
      DATA AE14CS( 2) /   -.0864811785 5259871E0 /
      DATA AE14CS( 3) /    .0072241015 4374659E0 /
      DATA AE14CS( 4) /   -.0008097559 4575573E0 /
      DATA AE14CS( 5) /    .0001099913 4432661E0 /
      DATA AE14CS( 6) /   -.0000171733 2998937E0 /
      DATA AE14CS( 7) /    .0000029856 2751447E0 /
      DATA AE14CS( 8) /   -.0000005659 6491457E0 /
      DATA AE14CS( 9) /    .0000001152 6808397E0 /
      DATA AE14CS(10) /   -.0000000249 5030440E0 /
      DATA AE14CS(11) /    .0000000056 9232420E0 /
      DATA AE14CS(12) /   -.0000000013 5995766E0 /
      DATA AE14CS(13) /    .0000000003 3846628E0 /
      DATA AE14CS(14) /   -.0000000000 8737853E0 /
      DATA AE14CS(15) /    .0000000000 2331588E0 /
      DATA AE14CS(16) /   -.0000000000 0641148E0 /
      DATA AE14CS(17) /    .0000000000 0181224E0 /
      DATA AE14CS(18) /   -.0000000000 0052538E0 /
      DATA AE14CS(19) /    .0000000000 0015592E0 /
      DATA AE14CS(20) /   -.0000000000 0004729E0 /
      DATA AE14CS(21) /    .0000000000 0001463E0 /
      DATA AE14CS(22) /   -.0000000000 0000461E0 /
      DATA AE14CS(23) /    .0000000000 0000148E0 /
      DATA AE14CS(24) /   -.0000000000 0000048E0 /
      DATA AE14CS(25) /    .0000000000 0000016E0 /
      DATA AE14CS(26) /   -.0000000000 0000005E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  E1
      IF (FIRST) THEN
         ETA = 0.1*R1MACH(3)
         NTAE11 = INITS (AE11CS, 39, ETA)
         NTAE12 = INITS (AE12CS, 25, ETA)
         NTE11 = INITS (E11CS, 19, ETA)
         NTE12 = INITS (E12CS, 16, ETA)
         NTAE13 = INITS (AE13CS, 25, ETA)
         NTAE14 = INITS (AE14CS, 26, ETA)
C
         XMAXT = -LOG (R1MACH(1))
         XMAX = XMAXT - LOG(XMAXT)
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GT.(-10.)) GO TO 20
C
C E1(X) = -EI(-X) FOR X .LE. -10.
C
      E1 = EXP(-X)/X * (1.+CSEVL (20./X+1., AE11CS, NTAE11))
      RETURN
C
 20   IF (X.GT.(-4.0)) GO TO 30
C
C E1(X) = -EI(-X) FOR -10. .LT. X .LE. -4.
C
      E1 = EXP(-X)/X * (1.+CSEVL ((40./X+7.)/3., AE12CS, NTAE12))
      RETURN
C
 30   IF (X.GT.(-1.0)) GO TO 40
C
C E1(X) = -EI(-X) FOR -4. .LT. X .LE. -1.
C
      E1 = -LOG(ABS(X)) + CSEVL ((2.*X+5.)/3., E11CS, NTE11)
      RETURN
C
 40   IF (X.GT.1.) GO TO 50
      IF (X .EQ. 0.) CALL XERMSG ('SLATEC', 'E1', 'X IS 0', 2, 2)
C
C E1(X) = -EI(-X) FOR -1. .LT. X .LE. 1.,  X .NE. 0.
C
      E1 = (-LOG(ABS(X)) - 0.6875 + X) + CSEVL (X, E12CS, NTE12)
      RETURN
C
 50   IF (X.GT.4.) GO TO 60
C
C E1(X) = -EI(-X) FOR 1. .LT. X .LE. 4.
C
      E1 = EXP(-X)/X * (1.+CSEVL ((8./X-5.)/3., AE13CS, NTAE13))
      RETURN
C
 60   IF (X.GT.XMAX) GO TO 70
C
C E1(X) = -EI(-X) FOR 4. .LT. X .LE. XMAX
C
      E1 = EXP(-X)/X * (1. + CSEVL (8./X-1., AE14CS, NTAE14))
      RETURN
C
C E1(X) = -EI(-X) FOR X .GT. XMAX
C
 70   CALL XERMSG ('SLATEC', 'E1', 'X SO BIG E1 UNDERFLOWS', 1, 1)
      E1 = 0.
      RETURN
C
      END
