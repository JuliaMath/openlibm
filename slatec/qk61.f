*DECK QK61
      SUBROUTINE QK61 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK61
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK61-S, DQK61-D)
C***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rule
C        Standard fortran subroutine
C        Real version
C
C
C        PARAMETERS
C         ON ENTRY
C           F      - Real
C                    Function subprogram defining the integrand
C                    function F(X). The actual name for F needs to be
C                    declared E X T E R N A L in the calling program.
C
C           A      - Real
C                    Lower limit of integration
C
C           B      - Real
C                    Upper limit of integration
C
C         ON RETURN
C           RESULT - Real
C                    Approximation to the integral I
C                    RESULT is computed by applying the 61-point
C                    Kronrod rule (RESK) obtained by optimal addition of
C                    abscissae to the 30-point Gauss rule (RESG).
C
C           ABSERR - Real
C                    Estimate of the modulus of the absolute error,
C                    which should equal or exceed ABS(I-RESULT)
C
C           RESABS - Real
C                    Approximation to the integral J
C
C           RESASC - Real
C                    Approximation to the integral of ABS(F-I/(B-A))
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK61
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(30),FV2(30),XGK(31),WGK(31),WG(15)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
C           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
C           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
C                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
C                   GAUSS RULE
C                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
C                   TO THE 30-POINT GAUSS RULE
C
C           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
C
C           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1   XGK(9),XGK(10)/
     2     0.9994844100504906E+00,     0.9968934840746495E+00,
     3     0.9916309968704046E+00,     0.9836681232797472E+00,
     4     0.9731163225011263E+00,     0.9600218649683075E+00,
     5     0.9443744447485600E+00,     0.9262000474292743E+00,
     6     0.9055733076999078E+00,     0.8825605357920527E+00/
      DATA XGK(11),XGK(12),XGK(13),XGK(14),XGK(15),XGK(16),
     1  XGK(17),XGK(18),XGK(19),XGK(20)/
     2     0.8572052335460611E+00,     0.8295657623827684E+00,
     3     0.7997278358218391E+00,     0.7677774321048262E+00,
     4     0.7337900624532268E+00,     0.6978504947933158E+00,
     5     0.6600610641266270E+00,     0.6205261829892429E+00,
     6     0.5793452358263617E+00,     0.5366241481420199E+00/
      DATA XGK(21),XGK(22),XGK(23),XGK(24),
     1  XGK(25),XGK(26),XGK(27),XGK(28),XGK(29),XGK(30),XGK(31)/
     2     0.4924804678617786E+00,     0.4470337695380892E+00,
     3     0.4004012548303944E+00,     0.3527047255308781E+00,
     4     0.3040732022736251E+00,     0.2546369261678898E+00,
     5     0.2045251166823099E+00,     0.1538699136085835E+00,
     6     0.1028069379667370E+00,     0.5147184255531770E-01,
     7     0.0E+00                   /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10)/
     2     0.1389013698677008E-02,     0.3890461127099884E-02,
     3     0.6630703915931292E-02,     0.9273279659517763E-02,
     4     0.1182301525349634E-01,     0.1436972950704580E-01,
     5     0.1692088918905327E-01,     0.1941414119394238E-01,
     6     0.2182803582160919E-01,     0.2419116207808060E-01/
      DATA WGK(11),WGK(12),WGK(13),WGK(14),WGK(15),WGK(16),
     1  WGK(17),WGK(18),WGK(19),WGK(20)/
     2     0.2650995488233310E-01,     0.2875404876504129E-01,
     3     0.3090725756238776E-01,     0.3298144705748373E-01,
     4     0.3497933802806002E-01,     0.3688236465182123E-01,
     5     0.3867894562472759E-01,     0.4037453895153596E-01,
     6     0.4196981021516425E-01,     0.4345253970135607E-01/
      DATA WGK(21),WGK(22),WGK(23),WGK(24),
     1  WGK(25),WGK(26),WGK(27),WGK(28),WGK(29),WGK(30),WGK(31)/
     2     0.4481480013316266E-01,     0.4605923827100699E-01,
     3     0.4718554656929915E-01,     0.4818586175708713E-01,
     4     0.4905543455502978E-01,     0.4979568342707421E-01,
     5     0.5040592140278235E-01,     0.5088179589874961E-01,
     6     0.5122154784925877E-01,     0.5142612853745903E-01,
     7     0.5149472942945157E-01/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.7968192496166606E-02,     0.1846646831109096E-01,
     2     0.2878470788332337E-01,     0.3879919256962705E-01,
     3     0.4840267283059405E-01,     0.5749315621761907E-01,
     4     0.6597422988218050E-01,     0.7375597473770521E-01/
      DATA WG(9),WG(10),WG(11),WG(12),WG(13),WG(14),WG(15)/
     1     0.8075589522942022E-01,     0.8689978720108298E-01,
     2     0.9212252223778613E-01,     0.9636873717464426E-01,
     3     0.9959342058679527E-01,     0.1017623897484055E+00,
     4     0.1028526528935588E+00/
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 30-POINT GAUSS RULE
C           RESK   - RESULT OF THE 61-POINT KRONROD RULE
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
C                    OVER (A,B), I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK61
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(B+A)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
C           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0E+00
      FC = F(CENTR)
      RESK = WGK(31)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,15
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J=1,15
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  15    CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(31)*ABS(FC-RESKH)
      DO 20 J=1,30
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.0E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
