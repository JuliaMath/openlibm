*DECK QK15I
      SUBROUTINE QK15I (F, BOUN, INF, A, B, RESULT, ABSERR, RESABS,
     +   RESASC)
C***BEGIN PROLOGUE  QK15I
C***PURPOSE  The original (infinite integration range is mapped
C            onto the interval (0,1) and (A,B) is a part of (0,1).
C            it is the purpose to compute
C            I = Integral of transformed integrand over (A,B),
C            J = Integral of ABS(Transformed Integrand) over (A,B).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A3A2, H2A4A2
C***TYPE      SINGLE PRECISION (QK15I-S, DQK15I-D)
C***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration Rule
C           Standard Fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              BOUN   - Real
C                       Finite bound of original integration
C                       Range (SET TO ZERO IF INF = +2)
C
C              INF    - Integer
C                       If INF = -1, the original interval is
C                                   (-INFINITY,BOUND),
C                       If INF = +1, the original interval is
C                                   (BOUND,+INFINITY),
C                       If INF = +2, the original interval is
C                                   (-INFINITY,+INFINITY) AND
C                       The integral is computed as the sum of two
C                       integrals, one over (-INFINITY,0) and one over
C                       (0,+INFINITY).
C
C              A      - Real
C                       Lower limit for integration over subrange
C                       of (0,1)
C
C              B      - Real
C                       Upper limit for integration over subrange
C                       of (0,1)
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       Result is computed by applying the 15-POINT
C                       KRONROD RULE(RESK) obtained by optimal addition
C                       of abscissae to the 7-POINT GAUSS RULE(RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the absolute error,
C                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of
C                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK15I
C
      REAL A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,
     1  DINF,R1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,
     2  FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,
     3  UFLOW,WG,WGK,XGK
      INTEGER INF,J
      EXTERNAL F
C
      DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
C           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
C           THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
C                    TO THE ABSCISSAE XGK(2), XGK(4), ...
C                    WG(1), WG(3), ... ARE SET TO ZERO.
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),
     1  XGK(8)/
     2     0.9914553711208126E+00,     0.9491079123427585E+00,
     3     0.8648644233597691E+00,     0.7415311855993944E+00,
     4     0.5860872354676911E+00,     0.4058451513773972E+00,
     5     0.2077849550078985E+00,     0.0000000000000000E+00/
C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),
     1  WGK(8)/
     2     0.2293532201052922E-01,     0.6309209262997855E-01,
     3     0.1047900103222502E+00,     0.1406532597155259E+00,
     4     0.1690047266392679E+00,     0.1903505780647854E+00,
     5     0.2044329400752989E+00,     0.2094821410847278E+00/
C
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.0000000000000000E+00,     0.1294849661688697E+00,
     2     0.0000000000000000E+00,     0.2797053914892767E+00,
     3     0.0000000000000000E+00,     0.3818300505051189E+00,
     4     0.0000000000000000E+00,     0.4179591836734694E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC*  - ABSCISSA
C           TABSC* - TRANSFORMED ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
C                    INTEGRAND OVER (A,B), I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK15I
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
      DINF = MIN(1,INF)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      TABSC1 = BOUN+DINF*(0.1E+01-CENTR)/CENTR
      FVAL1 = F(TABSC1)
      IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
      FC = (FVAL1/CENTR)/CENTR
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ERROR.
C
      RESG = WG(8)*FC
      RESK = WGK(8)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
        ABSC = HLGTH*XGK(J)
        ABSC1 = CENTR-ABSC
        ABSC2 = CENTR+ABSC
        TABSC1 = BOUN+DINF*(0.1E+01-ABSC1)/ABSC1
        TABSC2 = BOUN+DINF*(0.1E+01-ABSC2)/ABSC2
        FVAL1 = F(TABSC1)
        FVAL2 = F(TABSC2)
        IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
        IF(INF.EQ.2) FVAL2 = FVAL2+F(-TABSC2)
        FVAL1 = (FVAL1/ABSC1)/ABSC1
        FVAL2 = (FVAL2/ABSC2)/ABSC2
        FV1(J) = FVAL1
        FV2(J) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(J)*FSUM
        RESABS = RESABS+WGK(J)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESASC = RESASC*HLGTH
      RESABS = RESABS*HLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.E0) ABSERR = RESASC*
     1 MIN(0.1E+01,(0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1 ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
