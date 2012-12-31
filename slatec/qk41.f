*DECK QK41
      SUBROUTINE QK41 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK41
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK41-S, DQK41-D)
C***KEYWORDS  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
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
C                       RESULT is computed by applying the 41-POINT
C                       GAUSS-KRONROD RULE (RESK) obtained by optimal
C                       addition of abscissae to the 20-POINT GAUSS
C                       RULE (RESG).
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
C***END PROLOGUE  QK41
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,
     2  RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     3  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(20),FV2(20),XGK(21),WGK(21),WG(10)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 20-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE
C
C           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1  XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14),XGK(15),
     2  XGK(16),XGK(17),XGK(18),XGK(19),XGK(20),XGK(21)/
     3     0.9988590315882777E+00,   0.9931285991850949E+00,
     4     0.9815078774502503E+00,   0.9639719272779138E+00,
     5     0.9408226338317548E+00,   0.9122344282513259E+00,
     6     0.8782768112522820E+00,   0.8391169718222188E+00,
     7     0.7950414288375512E+00,   0.7463319064601508E+00,
     8     0.6932376563347514E+00,   0.6360536807265150E+00,
     9     0.5751404468197103E+00,   0.5108670019508271E+00,
     1     0.4435931752387251E+00,   0.3737060887154196E+00,
     2     0.3016278681149130E+00,   0.2277858511416451E+00,
     3     0.1526054652409227E+00,   0.7652652113349733E-01,
     4     0.0E+00               /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14),WGK(15),WGK(16),
     2  WGK(17),WGK(18),WGK(19),WGK(20),WGK(21)/
     3     0.3073583718520532E-02,   0.8600269855642942E-02,
     4     0.1462616925697125E-01,   0.2038837346126652E-01,
     5     0.2588213360495116E-01,   0.3128730677703280E-01,
     6     0.3660016975820080E-01,   0.4166887332797369E-01,
     7     0.4643482186749767E-01,   0.5094457392372869E-01,
     8     0.5519510534828599E-01,   0.5911140088063957E-01,
     9     0.6265323755478117E-01,   0.6583459713361842E-01,
     1     0.6864867292852162E-01,   0.7105442355344407E-01,
     2     0.7303069033278667E-01,   0.7458287540049919E-01,
     3     0.7570449768455667E-01,   0.7637786767208074E-01,
     4     0.7660071191799966E-01/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8),WG(9),WG(10)/
     1     0.1761400713915212E-01,    0.4060142980038694E-01,
     2     0.6267204833410906E-01,    0.8327674157670475E-01,
     3     0.1019301198172404E+00,    0.1181945319615184E+00,
     4     0.1316886384491766E+00,    0.1420961093183821E+00,
     5     0.1491729864726037E+00,    0.1527533871307259E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
C                    TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK41
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0E+00
      FC = F(CENTR)
      RESK = WGK(21)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,10
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
      DO 15 J = 1,10
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
      RESASC = WGK(21)*ABS(FC-RESKH)
      DO 20 J=1,20
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
