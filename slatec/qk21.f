*DECK QK21
      SUBROUTINE QK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK21
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK21-S, DQK21-D)
C***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
C                       Declared E X T E R N A L in the driver program.
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
C                       RESULT is computed by applying the 21-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 10-POINT GAUSS RULE (RESG).
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
C***END PROLOGUE  QK21
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,WG,WGK,
     2  XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 10-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),
     1  XGK(8),XGK(9),XGK(10),XGK(11)/
     2         0.9956571630258081E+00,     0.9739065285171717E+00,
     3     0.9301574913557082E+00,     0.8650633666889845E+00,
     4     0.7808177265864169E+00,     0.6794095682990244E+00,
     5     0.5627571346686047E+00,     0.4333953941292472E+00,
     6     0.2943928627014602E+00,     0.1488743389816312E+00,
     7     0.0000000000000000E+00/
C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),
     1  WGK(8),WGK(9),WGK(10),WGK(11)/
     2     0.1169463886737187E-01,     0.3255816230796473E-01,
     3     0.5475589657435200E-01,     0.7503967481091995E-01,
     4     0.9312545458369761E-01,     0.1093871588022976E+00,
     5     0.1234919762620659E+00,     0.1347092173114733E+00,
     6     0.1427759385770601E+00,     0.1477391049013385E+00,
     7     0.1494455540029169E+00/
C
      DATA WG(1),WG(2),WG(3),WG(4),WG(5)/
     1     0.6667134430868814E-01,     0.1494513491505806E+00,
     2     0.2190863625159820E+00,     0.2692667193099964E+00,
     3     0.2955242247147529E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK21
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0E+00
      FC = F(CENTR)
      RESK = WGK(11)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
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
      DO 15 J = 1,5
        JTWM1 = 2*J-1
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
      RESASC = WGK(11)*ABS(FC-RESKH)
      DO 20 J=1,10
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
