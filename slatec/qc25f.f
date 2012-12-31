*DECK QC25F
      SUBROUTINE QC25F (F, A, B, OMEGA, INTEGR, NRMOM, MAXP1, KSAVE,
     +   RESULT, ABSERR, NEVAL, RESABS, RESASC, MOMCOM, CHEBMO)
C***BEGIN PROLOGUE  QC25F
C***PURPOSE  To compute the integral I=Integral of F(X) over (A,B)
C            Where W(X) = COS(OMEGA*X) Or (WX)=SIN(OMEGA*X)
C            and to compute J=Integral of ABS(F) over (A,B). For small
C            value of OMEGA or small intervals (A,B) 15-point GAUSS-
C            KRONROD Rule used. Otherwise generalized CLENSHAW-CURTIS us
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A2
C***TYPE      SINGLE PRECISION (QC25F-S, DQC25F-D)
C***KEYWORDS  CLENSHAW-CURTIS METHOD, GAUSS-KRONROD RULES,
C             INTEGRATION RULES FOR FUNCTIONS WITH COS OR SIN FACTOR,
C             QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rules for functions with COS or SIN factor
C        Standard fortran subroutine
C        Real version
C
C        PARAMETERS
C         ON ENTRY
C           F      - Real
C                    Function subprogram defining the integrand
C                    function F(X). The actual name for F needs to
C                    be declared E X T E R N A L in the calling program.
C
C           A      - Real
C                    Lower limit of integration
C
C           B      - Real
C                    Upper limit of integration
C
C           OMEGA  - Real
C                    Parameter in the WEIGHT function
C
C           INTEGR - Integer
C                    Indicates which WEIGHT function is to be used
C                       INTEGR = 1   W(X) = COS(OMEGA*X)
C                       INTEGR = 2   W(X) = SIN(OMEGA*X)
C
C           NRMOM  - Integer
C                    The length of interval (A,B) is equal to the length
C                    of the original integration interval divided by
C                    2**NRMOM (we suppose that the routine is used in an
C                    adaptive integration process, otherwise set
C                    NRMOM = 0). NRMOM must be zero at the first call.
C
C           MAXP1  - Integer
C                    Gives an upper bound on the number of Chebyshev
C                    moments which can be stored, i.e. for the
C                    intervals of lengths ABS(BB-AA)*2**(-L),
C                    L = 0,1,2, ..., MAXP1-2.
C
C           KSAVE  - Integer
C                    Key which is one when the moments for the
C                    current interval have been computed
C
C         ON RETURN
C           RESULT - Real
C                    Approximation to the integral I
C
C           ABSERR - Real
C                    Estimate of the modulus of the absolute
C                    error, which should equal or exceed ABS(I-RESULT)
C
C           NEVAL  - Integer
C                    Number of integrand evaluations
C
C           RESABS - Real
C                    Approximation to the integral J
C
C           RESASC - Real
C                    Approximation to the integral of ABS(F-I/(B-A))
C
C         ON ENTRY AND RETURN
C           MOMCOM - Integer
C                    For each interval length we need to compute the
C                    Chebyshev moments. MOMCOM counts the number of
C                    intervals for which these moments have already been
C                    computed. If NRMOM.LT.MOMCOM or KSAVE = 1, the
C                    Chebyshev moments for the interval (A,B) have
C                    already been computed and stored, otherwise we
C                    compute them and we increase MOMCOM.
C
C           CHEBMO - Real
C                    Array of dimension at least (MAXP1,25) containing
C                    the modified Chebyshev moments for the first MOMCOM
C                    MOMCOM interval lengths
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QCHEB, QK15W, QWGTF, R1MACH, SGTSL
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QC25F
C
      REAL A,ABSERR,AC,AN,AN2,AS,ASAP,ASS,B,CENTR,CHEBMO,
     1  CHEB12,CHEB24,CONC,CONS,COSPAR,D,QWGTF,
     2  D1,R1MACH,D2,ESTC,ESTS,F,FVAL,HLGTH,OFLOW,OMEGA,PARINT,PAR2,
     3  PAR22,P2,P3,P4,RESABS,RESASC,RESC12,RESC24,RESS12,RESS24,
     4  RESULT,SINPAR,V,X
      INTEGER I,IERS,INTEGR,ISYM,J,K,KSAVE,M,MAXP1,MOMCOM,NEVAL,
     1  NOEQU,NOEQ1,NRMOM
C
      DIMENSION CHEBMO(MAXP1,25),CHEB12(13),CHEB24(25),D(25),D1(25),
     1  D2(25),FVAL(25),V(28),X(11)
C
      EXTERNAL F, QWGTF
C
C           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
C           K = 1, ...,11, TO BE USED FOR THE CHEBYSHEV EXPANSION OF F
C
      SAVE X
      DATA X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),
     1  X(10),X(11)/
     2     0.9914448613738104E+00,     0.9659258262890683E+00,
     3     0.9238795325112868E+00,     0.8660254037844386E+00,
     4     0.7933533402912352E+00,     0.7071067811865475E+00,
     5     0.6087614290087206E+00,     0.5000000000000000E+00,
     6     0.3826834323650898E+00,     0.2588190451025208E+00,
     7     0.1305261922200516E+00/
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
C           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
C                    (B-A)*0.5*COS(K*PI/12) + (B+A)*0.5,
C                    K = 0, ..., 24
C           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
C                    OF DEGREE 12, FOR THE FUNCTION F, IN THE
C                    INTERVAL (A,B)
C           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
C                    OF DEGREE 24, FOR THE FUNCTION F, IN THE
C                    INTERVAL (A,B)
C           RESC12 - APPROXIMATION TO THE INTEGRAL OF
C                    COS(0.5*(B-A)*OMEGA*X)*F(0.5*(B-A)*X+0.5*(B+A))
C                    OVER (-1,+1), USING THE CHEBYSHEV SERIES
C                    EXPANSION OF DEGREE 12
C           RESC24 - APPROXIMATION TO THE SAME INTEGRAL, USING THE
C                    CHEBYSHEV SERIES EXPANSION OF DEGREE 24
C           RESS12 - THE ANALOGUE OF RESC12 FOR THE SINE
C           RESS24 - THE ANALOGUE OF RESC24 FOR THE SINE
C
C
C           MACHINE DEPENDENT CONSTANT
C           --------------------------
C
C           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QC25F
      OFLOW = R1MACH(2)
C
      CENTR = 0.5E+00*(B+A)
      HLGTH = 0.5E+00*(B-A)
      PARINT = OMEGA*HLGTH
C
C           COMPUTE THE INTEGRAL USING THE 15-POINT GAUSS-KRONROD
C           FORMULA IF THE VALUE OF THE PARAMETER IN THE INTEGRAND
C           IS SMALL.
C
      IF(ABS(PARINT).GT.0.2E+01) GO TO 10
      CALL QK15W(F,QWGTF,OMEGA,P2,P3,P4,INTEGR,A,B,RESULT,
     1  ABSERR,RESABS,RESASC)
      NEVAL = 15
      GO TO 170
C
C           COMPUTE THE INTEGRAL USING THE GENERALIZED CLENSHAW-
C           CURTIS METHOD.
C
   10 CONC = HLGTH*COS(CENTR*OMEGA)
      CONS = HLGTH*SIN(CENTR*OMEGA)
      RESASC = OFLOW
      NEVAL = 25
C
C           CHECK WHETHER THE CHEBYSHEV MOMENTS FOR THIS INTERVAL
C           HAVE ALREADY BEEN COMPUTED.
C
      IF(NRMOM.LT.MOMCOM.OR.KSAVE.EQ.1) GO TO 120
C
C           COMPUTE A NEW SET OF CHEBYSHEV MOMENTS.
C
      M = MOMCOM+1
      PAR2 = PARINT*PARINT
      PAR22 = PAR2+0.2E+01
      SINPAR = SIN(PARINT)
      COSPAR = COS(PARINT)
C
C           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO COSINE.
C
      V(1) = 0.2E+01*SINPAR/PARINT
      V(2) = (0.8E+01*COSPAR+(PAR2+PAR2-0.8E+01)*SINPAR/
     1  PARINT)/PAR2
      V(3) = (0.32E+02*(PAR2-0.12E+02)*COSPAR+(0.2E+01*
     1  ((PAR2-0.80E+02)*PAR2+0.192E+03)*SINPAR)/
     2  PARINT)/(PAR2*PAR2)
      AC = 0.8E+01*COSPAR
      AS = 0.24E+02*PARINT*SINPAR
      IF(ABS(PARINT).GT.0.24E+02) GO TO 30
C
C           COMPUTE THE CHEBYSHEV MOMENTS AS THE
C           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1
C           INITIAL VALUE (V(3)) AND 1 END VALUE (COMPUTED
C           USING AN ASYMPTOTIC FORMULA).
C
      NOEQU = 25
      NOEQ1 = NOEQU-1
      AN = 0.6E+01
      DO 20 K = 1,NOEQ1
        AN2 = AN*AN
        D(K) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
        D2(K) = (AN-0.1E+01)*(AN-0.2E+01)*PAR2
        D1(K+1) = (AN+0.3E+01)*(AN+0.4E+01)*PAR2
        V(K+3) = AS-(AN2-0.4E+01)*AC
        AN = AN+0.2E+01
   20 CONTINUE
      AN2 = AN*AN
      D(NOEQU) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
      V(NOEQU+3) = AS-(AN2-0.4E+01)*AC
      V(4) = V(4)-0.56E+02*PAR2*V(3)
      ASS = PARINT*SINPAR
      ASAP = (((((0.210E+03*PAR2-0.1E+01)*COSPAR-(0.105E+03*PAR2
     1  -0.63E+02)*ASS)/AN2-(0.1E+01-0.15E+02*PAR2)*COSPAR
     2  +0.15E+02*ASS)/AN2-COSPAR+0.3E+01*ASS)/AN2-COSPAR)/AN2
      V(NOEQU+3) = V(NOEQU+3)-0.2E+01*ASAP*PAR2*(AN-0.1E+01)*
     1   (AN-0.2E+01)
C
C           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN
C           ELIMINATION WITH PARTIAL PIVOTING.
C
      CALL SGTSL(NOEQU,D1,D,D2,V(4),IERS)
      GO TO 50
C
C           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD
C           RECURSION.
C
   30 AN = 0.4E+01
      DO 40 I = 4,13
        AN2 = AN*AN
        V(I) = ((AN2-0.4E+01)*(0.2E+01*(PAR22-AN2-AN2)*V(I-1)-AC)
     1  +AS-PAR2*(AN+0.1E+01)*(AN+0.2E+01)*V(I-2))/
     2  (PAR2*(AN-0.1E+01)*(AN-0.2E+01))
        AN = AN+0.2E+01
   40 CONTINUE
   50 DO 60 J = 1,13
        CHEBMO(M,2*J-1) = V(J)
   60 CONTINUE
C
C           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO SINE.
C
      V(1) = 0.2E+01*(SINPAR-PARINT*COSPAR)/PAR2
      V(2) = (0.18E+02-0.48E+02/PAR2)*SINPAR/PAR2
     1  +(-0.2E+01+0.48E+02/PAR2)*COSPAR/PARINT
      AC = -0.24E+02*PARINT*COSPAR
      AS = -0.8E+01*SINPAR
      IF(ABS(PARINT).GT.0.24E+02) GO TO 80
C
C           COMPUTE THE CHEBYSHEV MOMENTS AS THE
C           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1
C           INITIAL VALUE (V(2)) AND 1 END VALUE (COMPUTED
C           USING AN ASYMPTOTIC FORMULA).
C
      AN = 0.5E+01
      DO 70 K = 1,NOEQ1
        AN2 = AN*AN
        D(K) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
        D2(K) = (AN-0.1E+01)*(AN-0.2E+01)*PAR2
        D1(K+1) = (AN+0.3E+01)*(AN+0.4E+01)*PAR2
        V(K+2) = AC+(AN2-0.4E+01)*AS
        AN = AN+0.2E+01
   70 CONTINUE
      AN2 = AN*AN
      D(NOEQU) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
      V(NOEQU+2) = AC+(AN2-0.4E+01)*AS
      V(3) = V(3)-0.42E+02*PAR2*V(2)
      ASS = PARINT*COSPAR
      ASAP = (((((0.105E+03*PAR2-0.63E+02)*ASS+(0.210E+03*PAR2
     1  -0.1E+01)*SINPAR)/AN2+(0.15E+02*PAR2-0.1E+01)*SINPAR-
     2  0.15E+02*ASS)/AN2-0.3E+01*ASS-SINPAR)/AN2-SINPAR)/AN2
      V(NOEQU+2) = V(NOEQU+2)-0.2E+01*ASAP*PAR2*(AN-0.1E+01)
     1  *(AN-0.2E+01)
C
C           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN
C           ELIMINATION WITH PARTIAL PIVOTING.
C
      CALL SGTSL(NOEQU,D1,D,D2,V(3),IERS)
      GO TO 100
C
C           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF
C           FORWARD RECURSION.
C
   80 AN = 0.3E+01
      DO 90 I = 3,12
        AN2 = AN*AN
        V(I) = ((AN2-0.4E+01)*(0.2E+01*(PAR22-AN2-AN2)*V(I-1)+AS)
     1  +AC-PAR2*(AN+0.1E+01)*(AN+0.2E+01)*V(I-2))
     2  /(PAR2*(AN-0.1E+01)*(AN-0.2E+01))
        AN = AN+0.2E+01
   90 CONTINUE
  100 DO 110 J = 1,12
        CHEBMO(M,2*J) = V(J)
  110 CONTINUE
  120 IF (NRMOM.LT.MOMCOM) M = NRMOM+1
       IF (MOMCOM.LT.MAXP1-1.AND.NRMOM.GE.MOMCOM) MOMCOM = MOMCOM+1
C
C           COMPUTE THE COEFFICIENTS OF THE CHEBYSHEV EXPANSIONS
C           OF DEGREES 12 AND 24 OF THE FUNCTION F.
C
      FVAL(1) = 0.5E+00*F(CENTR+HLGTH)
      FVAL(13) = F(CENTR)
      FVAL(25) = 0.5E+00*F(CENTR-HLGTH)
      DO 130 I = 2,12
        ISYM = 26-I
        FVAL(I) = F(HLGTH*X(I-1)+CENTR)
        FVAL(ISYM) = F(CENTR-HLGTH*X(I-1))
  130 CONTINUE
      CALL QCHEB(X,FVAL,CHEB12,CHEB24)
C
C           COMPUTE THE INTEGRAL AND ERROR ESTIMATES.
C
      RESC12 = CHEB12(13)*CHEBMO(M,13)
      RESS12 = 0.0E+00
      K = 11
      DO 140 J = 1,6
        RESC12 = RESC12+CHEB12(K)*CHEBMO(M,K)
        RESS12 = RESS12+CHEB12(K+1)*CHEBMO(M,K+1)
        K = K-2
  140 CONTINUE
      RESC24 = CHEB24(25)*CHEBMO(M,25)
      RESS24 = 0.0E+00
      RESABS = ABS(CHEB24(25))
      K = 23
      DO 150 J = 1,12
        RESC24 = RESC24+CHEB24(K)*CHEBMO(M,K)
        RESS24 = RESS24+CHEB24(K+1)*CHEBMO(M,K+1)
        RESABS = ABS(CHEB24(K))+ABS(CHEB24(K+1))
        K = K-2
  150 CONTINUE
      ESTC = ABS(RESC24-RESC12)
      ESTS = ABS(RESS24-RESS12)
      RESABS = RESABS*ABS(HLGTH)
      IF(INTEGR.EQ.2) GO TO 160
      RESULT = CONC*RESC24-CONS*RESS24
      ABSERR = ABS(CONC*ESTC)+ABS(CONS*ESTS)
      GO TO 170
  160 RESULT = CONC*RESS24+CONS*RESC24
      ABSERR = ABS(CONC*ESTS)+ABS(CONS*ESTC)
  170 RETURN
      END
