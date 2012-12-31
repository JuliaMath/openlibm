*DECK QK51
      SUBROUTINE QK51 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK51
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK51-S, DQK51-D)
C***KEYWORDS  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subroutine defining the integrand
C                       function F(X). The actual name for F needs to be
C                       declared E X T E R N A L in the calling program.
C
C              A      - Real
C                       Lower limit of integration
C
C              B      - Real
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       RESULT is computed by applying the 51-point
C                       Kronrod rule (RESK) obtained by optimal addition
C                       of abscissae to the 25-point Gauss rule (RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK51
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(25),FV2(25),XGK(26),WGK(26),WG(13)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 25-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1  XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14)/
     2     0.9992621049926098E+00,   0.9955569697904981E+00,
     3     0.9880357945340772E+00,   0.9766639214595175E+00,
     4     0.9616149864258425E+00,   0.9429745712289743E+00,
     5     0.9207471152817016E+00,   0.8949919978782754E+00,
     6     0.8658470652932756E+00,   0.8334426287608340E+00,
     7     0.7978737979985001E+00,   0.7592592630373576E+00,
     8     0.7177664068130844E+00,   0.6735663684734684E+00/
       DATA XGK(15),XGK(16),XGK(17),XGK(18),XGK(19),XGK(20),XGK(21),
     1  XGK(22),XGK(23),XGK(24),XGK(25),XGK(26)/
     2     0.6268100990103174E+00,   0.5776629302412230E+00,
     3     0.5263252843347192E+00,   0.4730027314457150E+00,
     4     0.4178853821930377E+00,   0.3611723058093878E+00,
     5     0.3030895389311078E+00,   0.2438668837209884E+00,
     6     0.1837189394210489E+00,   0.1228646926107104E+00,
     7     0.6154448300568508E-01,   0.0E+00               /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14)/
     2     0.1987383892330316E-02,   0.5561932135356714E-02,
     3     0.9473973386174152E-02,   0.1323622919557167E-01,
     4     0.1684781770912830E-01,   0.2043537114588284E-01,
     5     0.2400994560695322E-01,   0.2747531758785174E-01,
     6     0.3079230016738749E-01,   0.3400213027432934E-01,
     7     0.3711627148341554E-01,   0.4008382550403238E-01,
     8     0.4287284502017005E-01,   0.4550291304992179E-01/
       DATA WGK(15),WGK(16),WGK(17),WGK(18),WGK(19),WGK(20),WGK(21)
     1  ,WGK(22),WGK(23),WGK(24),WGK(25),WGK(26)/
     2     0.4798253713883671E-01,   0.5027767908071567E-01,
     3     0.5236288580640748E-01,   0.5425112988854549E-01,
     4     0.5595081122041232E-01,   0.5743711636156783E-01,
     5     0.5868968002239421E-01,   0.5972034032417406E-01,
     6     0.6053945537604586E-01,   0.6112850971705305E-01,
     7     0.6147118987142532E-01,   0.6158081806783294E-01/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8),WG(9),
     1  WG(10),WG(11),WG(12),WG(13)/
     2     0.1139379850102629E-01,    0.2635498661503214E-01,
     3     0.4093915670130631E-01,    0.5490469597583519E-01,
     4     0.6803833381235692E-01,    0.8014070033500102E-01,
     5     0.9102826198296365E-01,    0.1005359490670506E+00,
     6     0.1085196244742637E+00,    0.1148582591457116E+00,
     7     0.1194557635357848E+00,    0.1222424429903100E+00,
     8     0.1231760537267155E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK51
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(13)*FC
      RESK = WGK(26)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,12
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
      DO 15 J = 1,13
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(26)*ABS(FC-RESKH)
      DO 20 J=1,25
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
