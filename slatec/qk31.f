*DECK QK31
      SUBROUTINE QK31 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK31
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK31-S, DQK31-D)
C***KEYWORDS  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
C                       Declared E X T E R N A L in the calling program.
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
C                       RESULT is computed by applying the 31-POINT
C                       GAUSS-KRONROD RULE (RESK), obtained by optimal
C                       addition of abscissae to the 15-POINT GAUSS
C                       RULE (RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the modulus,
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
C***END PROLOGUE  QK31
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(15),FV2(15),XGK(16),WGK(16),WG(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 15-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1  XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14),XGK(15),
     2  XGK(16)/
     3     0.9980022986933971E+00,   0.9879925180204854E+00,
     4     0.9677390756791391E+00,   0.9372733924007059E+00,
     5     0.8972645323440819E+00,   0.8482065834104272E+00,
     6     0.7904185014424659E+00,   0.7244177313601700E+00,
     7     0.6509967412974170E+00,   0.5709721726085388E+00,
     8     0.4850818636402397E+00,   0.3941513470775634E+00,
     9     0.2991800071531688E+00,   0.2011940939974345E+00,
     1     0.1011420669187175E+00,   0.0E+00               /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14),WGK(15),
     2  WGK(16)/
     3     0.5377479872923349E-02,   0.1500794732931612E-01,
     4     0.2546084732671532E-01,   0.3534636079137585E-01,
     5     0.4458975132476488E-01,   0.5348152469092809E-01,
     6     0.6200956780067064E-01,   0.6985412131872826E-01,
     7     0.7684968075772038E-01,   0.8308050282313302E-01,
     8     0.8856444305621177E-01,   0.9312659817082532E-01,
     9     0.9664272698362368E-01,   0.9917359872179196E-01,
     1     0.1007698455238756E+00,   0.1013300070147915E+00/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.3075324199611727E-01,   0.7036604748810812E-01,
     2     0.1071592204671719E+00,   0.1395706779261543E+00,
     3     0.1662692058169939E+00,   0.1861610000155622E+00,
     4     0.1984314853271116E+00,   0.2025782419255613E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK31
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(8)*FC
      RESK = WGK(16)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
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
      DO 15 J = 1,8
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
      RESASC = WGK(16)*ABS(FC-RESKH)
      DO 20 J=1,15
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
