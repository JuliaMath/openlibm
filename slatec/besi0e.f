*DECK BESI0E
      FUNCTION BESI0E (X)
C***BEGIN PROLOGUE  BESI0E
C***PURPOSE  Compute the exponentially scaled modified (hyperbolic)
C            Bessel function of the first kind of order zero.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10B1
C***TYPE      SINGLE PRECISION (BESI0E-S, DBSI0E-D)
C***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
C             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
C             ORDER ZERO, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BESI0E(X) calculates the exponentially scaled modified (hyperbolic)
C Bessel function of the first kind of order zero for real argument X;
C i.e., EXP(-ABS(X))*I0(X).
C
C
C Series for BI0        on the interval  0.          to  9.00000D+00
C                                        with weighted error   2.46E-18
C                                         log weighted error  17.61
C                               significant figures required  17.90
C                                    decimal places required  18.15
C
C
C Series for AI0        on the interval  1.25000D-01 to  3.33333D-01
C                                        with weighted error   7.87E-17
C                                         log weighted error  16.10
C                               significant figures required  14.69
C                                    decimal places required  16.76
C
C
C Series for AI02       on the interval  0.          to  1.25000D-01
C                                        with weighted error   3.79E-17
C                                         log weighted error  16.42
C                               significant figures required  14.86
C                                    decimal places required  17.09
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CSEVL, INITS, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890313  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  BESI0E
      DIMENSION BI0CS(12), AI0CS(21), AI02CS(22)
      LOGICAL FIRST
      SAVE BI0CS, AI0CS, AI02CS, NTI0, NTAI0, NTAI02, XSML, FIRST
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
      DATA AI0CS( 1) /    .0757599449 4023796E0 /
      DATA AI0CS( 2) /    .0075913808 1082334E0 /
      DATA AI0CS( 3) /    .0004153131 3389237E0 /
      DATA AI0CS( 4) /    .0000107007 6463439E0 /
      DATA AI0CS( 5) /   -.0000079011 7997921E0 /
      DATA AI0CS( 6) /   -.0000007826 1435014E0 /
      DATA AI0CS( 7) /    .0000002783 8499429E0 /
      DATA AI0CS( 8) /    .0000000082 5247260E0 /
      DATA AI0CS( 9) /   -.0000000120 4463945E0 /
      DATA AI0CS(10) /    .0000000015 5964859E0 /
      DATA AI0CS(11) /    .0000000002 2925563E0 /
      DATA AI0CS(12) /   -.0000000001 1916228E0 /
      DATA AI0CS(13) /    .0000000000 1757854E0 /
      DATA AI0CS(14) /    .0000000000 0112822E0 /
      DATA AI0CS(15) /   -.0000000000 0114684E0 /
      DATA AI0CS(16) /    .0000000000 0027155E0 /
      DATA AI0CS(17) /   -.0000000000 0002415E0 /
      DATA AI0CS(18) /   -.0000000000 0000608E0 /
      DATA AI0CS(19) /    .0000000000 0000314E0 /
      DATA AI0CS(20) /   -.0000000000 0000071E0 /
      DATA AI0CS(21) /    .0000000000 0000007E0 /
      DATA AI02CS( 1) /    .0544904110 1410882E0 /
      DATA AI02CS( 2) /    .0033691164 7825569E0 /
      DATA AI02CS( 3) /    .0000688975 8346918E0 /
      DATA AI02CS( 4) /    .0000028913 7052082E0 /
      DATA AI02CS( 5) /    .0000002048 9185893E0 /
      DATA AI02CS( 6) /    .0000000226 6668991E0 /
      DATA AI02CS( 7) /    .0000000033 9623203E0 /
      DATA AI02CS( 8) /    .0000000004 9406022E0 /
      DATA AI02CS( 9) /    .0000000000 1188914E0 /
      DATA AI02CS(10) /   -.0000000000 3149915E0 /
      DATA AI02CS(11) /   -.0000000000 1321580E0 /
      DATA AI02CS(12) /   -.0000000000 0179419E0 /
      DATA AI02CS(13) /    .0000000000 0071801E0 /
      DATA AI02CS(14) /    .0000000000 0038529E0 /
      DATA AI02CS(15) /    .0000000000 0001539E0 /
      DATA AI02CS(16) /   -.0000000000 0004151E0 /
      DATA AI02CS(17) /   -.0000000000 0000954E0 /
      DATA AI02CS(18) /    .0000000000 0000382E0 /
      DATA AI02CS(19) /    .0000000000 0000176E0 /
      DATA AI02CS(20) /   -.0000000000 0000034E0 /
      DATA AI02CS(21) /   -.0000000000 0000027E0 /
      DATA AI02CS(22) /    .0000000000 0000003E0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BESI0E
      IF (FIRST) THEN
         NTI0 = INITS (BI0CS, 12, 0.1*R1MACH(3))
         NTAI0 = INITS (AI0CS, 21, 0.1*R1MACH(3))
         NTAI02 = INITS (AI02CS, 22, 0.1*R1MACH(3))
         XSML = SQRT (4.5*R1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.3.0) GO TO 20
C
      BESI0E = 1.0 - X
      IF (Y.GT.XSML) BESI0E = EXP(-Y) * ( 2.75 +
     1  CSEVL (Y*Y/4.5-1.0, BI0CS, NTI0) )
      RETURN
C
 20   IF (Y.LE.8.) BESI0E = (.375 + CSEVL ((48./Y-11.)/5., AI0CS, NTAI0)
     1  ) / SQRT(Y)
      IF (Y.GT.8.) BESI0E = (.375 + CSEVL (16./Y-1., AI02CS, NTAI02))
     1  / SQRT(Y)
C
      RETURN
      END
