*DECK QC25C
      SUBROUTINE QC25C (F, A, B, C, RESULT, ABSERR, KRUL, NEVAL)
C***BEGIN PROLOGUE  QC25C
C***PURPOSE  To compute I = Integral of F*W over (A,B) with
C            error estimate, where W(X) = 1/(X-C)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A2, J4
C***TYPE      SINGLE PRECISION (QC25C-S, DQC25C-D)
C***KEYWORDS  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rules for the computation of CAUCHY
C        PRINCIPAL VALUE integrals
C        Standard fortran subroutine
C        Real version
C
C        PARAMETERS
C           F      - Real
C                    Function subprogram defining the integrand function
C                    F(X). The actual name for F needs to be declared
C                    E X T E R N A L  in the driver program.
C
C           A      - Real
C                    Left end point of the integration interval
C
C           B      - Real
C                    Right end point of the integration interval, B.GT.A
C
C           C      - Real
C                    Parameter in the WEIGHT function
C
C           RESULT - Real
C                    Approximation to the integral
C                    result is computed by using a generalized
C                    Clenshaw-Curtis method if C lies within ten percent
C                    of the integration interval. In the other case the
C                    15-point Kronrod rule obtained by optimal addition
C                    of abscissae to the 7-point Gauss rule, is applied.
C
C           ABSERR - Real
C                    Estimate of the modulus of the absolute error,
C                    which should equal or exceed ABS(I-RESULT)
C
C           KRUL   - Integer
C                    Key which is decreased by 1 if the 15-point
C                    Gauss-Kronrod scheme has been used
C
C           NEVAL  - Integer
C                    Number of integrand evaluations
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QCHEB, QK15W, QWGTC
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QC25C
C
      REAL A,ABSERR,AK22,AMOM0,AMOM1,AMOM2,B,C,CC,
     1  CENTR,CHEB12,CHEB24,QWGTC,F,FVAL,HLGTH,P2,P3,P4,
     2  RESABS,RESASC,RESULT,RES12,RES24,U,X
      INTEGER I,ISYM,K,KP,KRUL,NEVAL
C
      DIMENSION X(11),FVAL(25),CHEB12(13),CHEB24(25)
C
      EXTERNAL F, QWGTC
C
C           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24),
C           K = 1, ..., 11, TO BE USED FOR THE CHEBYSHEV SERIES
C           EXPANSION OF F
C
      SAVE X
      DATA X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10),
     1  X(11)/
     2     0.9914448613738104E+00,     0.9659258262890683E+00,
     3     0.9238795325112868E+00,     0.8660254037844386E+00,
     4     0.7933533402912352E+00,     0.7071067811865475E+00,
     5     0.6087614290087206E+00,     0.5000000000000000E+00,
     6     0.3826834323650898E+00,     0.2588190451025208E+00,
     7     0.1305261922200516E+00/
C
C           LIST OF MAJOR VARIABLES
C           ----------------------
C           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
C                    COS(K*PI/24),  K = 0, ..., 24
C           CHEB12 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
C                    FOR THE FUNCTION F, OF DEGREE 12
C           CHEB24 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
C                    FOR THE FUNCTION F, OF DEGREE 24
C           RES12  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
C                    TO THE USE OF CHEB12
C           RES24  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
C                    TO THE USE OF CHEB24
C           QWGTC - EXTERNAL FUNCTION SUBPROGRAM DEFINING
C                    THE WEIGHT FUNCTION
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           CENTR  - MID POINT OF THE INTERVAL
C
C
C           CHECK THE POSITION OF C.
C
C***FIRST EXECUTABLE STATEMENT  QC25C
      CC = (0.2E+01*C-B-A)/(B-A)
      IF(ABS(CC).LT.0.11E+01) GO TO 10
C
C           APPLY THE 15-POINT GAUSS-KRONROD SCHEME.
C
      KRUL = KRUL-1
      CALL QK15W(F,QWGTC,C,P2,P3,P4,KP,A,B,RESULT,ABSERR,
     1  RESABS,RESASC)
      NEVAL = 15
      IF (RESASC.EQ.ABSERR) KRUL = KRUL+1
      GO TO 50
C
C           USE THE GENERALIZED CLENSHAW-CURTIS METHOD.
C
   10 HLGTH = 0.5E+00*(B-A)
      CENTR = 0.5E+00*(B+A)
      NEVAL = 25
      FVAL(1) = 0.5E+00*F(HLGTH+CENTR)
      FVAL(13) = F(CENTR)
      FVAL(25) = 0.5E+00*F(CENTR-HLGTH)
      DO 20 I=2,12
        U = HLGTH*X(I-1)
        ISYM = 26-I
        FVAL(I) = F(U+CENTR)
        FVAL(ISYM) = F(CENTR-U)
   20 CONTINUE
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION.
C
      CALL QCHEB(X,FVAL,CHEB12,CHEB24)
C
C           THE MODIFIED CHEBYSHEV MOMENTS ARE COMPUTED
C           BY FORWARD RECURSION, USING AMOM0 AND AMOM1
C           AS STARTING VALUES.
C
      AMOM0 = LOG(ABS((0.1E+01-CC)/(0.1E+01+CC)))
      AMOM1 = 0.2E+01+CC*AMOM0
      RES12 = CHEB12(1)*AMOM0+CHEB12(2)*AMOM1
      RES24 = CHEB24(1)*AMOM0+CHEB24(2)*AMOM1
      DO 30 K=3,13
        AMOM2 = 0.2E+01*CC*AMOM1-AMOM0
        AK22 = (K-2)*(K-2)
        IF((K/2)*2.EQ.K) AMOM2 = AMOM2-0.4E+01/(AK22-0.1E+01)
        RES12 = RES12+CHEB12(K)*AMOM2
        RES24 = RES24+CHEB24(K)*AMOM2
        AMOM0 = AMOM1
        AMOM1 = AMOM2
   30 CONTINUE
      DO 40 K=14,25
        AMOM2 = 0.2E+01*CC*AMOM1-AMOM0
        AK22 = (K-2)*(K-2)
        IF((K/2)*2.EQ.K) AMOM2 = AMOM2-0.4E+01/
     1  (AK22-0.1E+01)
        RES24 = RES24+CHEB24(K)*AMOM2
        AMOM0 = AMOM1
        AMOM1 = AMOM2
   40 CONTINUE
      RESULT = RES24
      ABSERR = ABS(RES24-RES12)
   50 RETURN
      END
