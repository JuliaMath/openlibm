*DECK DQK15W
      SUBROUTINE DQK15W (F, W, P1, P2, P3, P4, KP, A, B, RESULT, ABSERR,
     +   RESABS, RESASC)
C***BEGIN PROLOGUE  DQK15W
C***PURPOSE  To compute I = Integral of F*W over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F*W) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A2
C***TYPE      DOUBLE PRECISION (QK15W-S, DQK15W-D)
C***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
C           Double precision version
C
C           PARAMETERS
C             ON ENTRY
C              F      - Double precision
C                       Function subprogram defining the integrand
C                       function F(X). The actual name for F needs to be
C                       declared E X T E R N A L in the driver program.
C
C              W      - Double precision
C                       Function subprogram defining the integrand
C                       WEIGHT function W(X). The actual name for W
C                       needs to be declared E X T E R N A L in the
C                       calling program.
C
C              P1, P2, P3, P4 - Double precision
C                       Parameters in the WEIGHT function
C
C              KP     - Integer
C                       Key for indicating the type of WEIGHT function
C
C              A      - Double precision
C                       Lower limit of integration
C
C              B      - Double precision
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Double precision
C                       Approximation to the integral I
C                       RESULT is computed by applying the 15-point
C                       Kronrod rule (RESK) obtained by optimal addition
C                       of abscissae to the 7-point Gauss rule (RESG).
C
C              ABSERR - Double precision
C                       Estimate of the modulus of the absolute error,
C                       which should equal or exceed ABS(I-RESULT)
C
C              RESABS - Double precision
C                       Approximation to the integral of ABS(F)
C
C              RESASC - Double precision
C                       Approximation to the integral of ABS(F-I/(B-A))
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  DQK15W
C
      DOUBLE PRECISION A,ABSC,ABSC1,ABSC2,ABSERR,B,CENTR,DHLGTH,
     1  D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,
     2  P1,P2,P3,P4,RESABS,RESASC,RESG,RESK,RESKH,RESULT,UFLOW,W,WG,WGK,
     3  XGK
      INTEGER J,JTW,JTWM1,KP
      EXTERNAL F, W
C
      DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(4)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT GAUSS-KRONROD RULE
C                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ... ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT GAUSS-KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/
     1     0.9914553711208126D+00,     0.9491079123427585D+00,
     2     0.8648644233597691D+00,     0.7415311855993944D+00,
     3     0.5860872354676911D+00,     0.4058451513773972D+00,
     4     0.2077849550078985D+00,     0.0000000000000000D+00/
C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/
     1     0.2293532201052922D-01,     0.6309209262997855D-01,
     2     0.1047900103222502D+00,     0.1406532597155259D+00,
     3     0.1690047266392679D+00,     0.1903505780647854D+00,
     4     0.2044329400752989D+00,     0.2094821410847278D+00/
C
      DATA WG(1),WG(2),WG(3),WG(4)/
     1     0.1294849661688697D+00,    0.2797053914892767D+00,
     2     0.3818300505051889D+00,    0.4179591836734694D+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC*  - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F*W OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  DQK15W
      EPMACH = D1MACH(4)
      UFLOW = D1MACH(1)
C
      CENTR = 0.5D+00*(A+B)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE
C           INTEGRAL, AND ESTIMATE THE ERROR.
C
      FC = F(CENTR)*W(CENTR,P1,P2,P3,P4,KP)
      RESG = WG(4)*FC
      RESK = WGK(8)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        ABSC1 = CENTR-ABSC
        ABSC2 = CENTR+ABSC
        FVAL1 = F(ABSC1)*W(ABSC1,P1,P2,P3,P4,KP)
        FVAL2 = F(ABSC2)*W(ABSC2,P1,P2,P3,P4,KP)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J=1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        ABSC1 = CENTR-ABSC
        ABSC2 = CENTR+ABSC
        FVAL1 = F(ABSC1)*W(ABSC1,P1,P2,P3,P4,KP)
        FVAL2 = F(ABSC2)*W(ABSC2,P1,P2,P3,P4,KP)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX((EPMACH*
     1  0.5D+02)*RESABS,ABSERR)
      RETURN
      END
