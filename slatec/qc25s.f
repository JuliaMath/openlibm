*DECK QC25S
      SUBROUTINE QC25S (F, A, B, BL, BR, ALFA, BETA, RI, RJ, RG, RH,
     +   RESULT, ABSERR, RESASC, INTEGR, NEV)
C***BEGIN PROLOGUE  QC25S
C***PURPOSE  To compute I = Integral of F*W over (BL,BR), with error
C            estimate, where the weight function W has a singular
C            behaviour of ALGEBRAICO-LOGARITHMIC type at the points
C            A and/or B. (BL,BR) is a part of (A,B).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A2A2
C***TYPE      SINGLE PRECISION (QC25S-S, DQC25S-D)
C***KEYWORDS  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rules for integrands having ALGEBRAICO-LOGARITHMIC
C        end point singularities
C        Standard fortran subroutine
C        Real version
C
C        PARAMETERS
C           F      - Real
C                    Function subprogram defining the integrand
C                    F(X). The actual name for F needs to be declared
C                    E X T E R N A L  in the driver program.
C
C           A      - Real
C                    Left end point of the original interval
C
C           B      - Real
C                    Right end point of the original interval, B.GT.A
C
C           BL     - Real
C                    Lower limit of integration, BL.GE.A
C
C           BR     - Real
C                    Upper limit of integration, BR.LE.B
C
C           ALFA   - Real
C                    PARAMETER IN THE WEIGHT FUNCTION
C
C           BETA   - Real
C                    Parameter in the weight function
C
C           RI,RJ,RG,RH - Real
C                    Modified CHEBYSHEV moments for the application
C                    of the generalized CLENSHAW-CURTIS
C                    method (computed in subroutine DQMOMO)
C
C           RESULT - Real
C                    Approximation to the integral
C                    RESULT is computed by using a generalized
C                    CLENSHAW-CURTIS method if B1 = A or BR = B.
C                    in all other cases the 15-POINT KRONROD
C                    RULE is applied, obtained by optimal addition of
C                    Abscissae to the 7-POINT GAUSS RULE.
C
C           ABSERR - Real
C                    Estimate of the modulus of the absolute error,
C                    which should equal or exceed ABS(I-RESULT)
C
C           RESASC - Real
C                    Approximation to the integral of ABS(F*W-I/(B-A))
C
C           INTEGR - Integer
C                    Which determines the weight function
C                    = 1   W(X) = (X-A)**ALFA*(B-X)**BETA
C                    = 2   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
C                    = 3   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
C                    = 4   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*
C                                 LOG(B-X)
C
C           NEV    - Integer
C                    Number of integrand evaluations
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QCHEB, QK15W, QWGTS
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QC25S
C
      REAL A,ABSERR,ALFA,B,BETA,BL,BR,CENTR,CHEB12,CHEB24,
     1  DC,F,FACTOR,FIX,FVAL,HLGTH,RESABS,RESASC,
     2  RESULT,RES12,RES24,RG,RH,RI,RJ,U,QWGTS,X
      INTEGER I,INTEGR,ISYM,NEV
C
      DIMENSION CHEB12(13),CHEB24(25),FVAL(25),RG(25),RH(25),RI(25),
     1  RJ(25),X(11)
C
      EXTERNAL F, QWGTS
C
C           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
C           K = 1, ..., 11, TO BE USED FOR THE COMPUTATION OF THE
C           CHEBYSHEV SERIES EXPANSION OF F.
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
C           -----------------------
C
C           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
C                    (BR-BL)*0.5*COS(K*PI/24)+(BR+BL)*0.5
C                    K = 0, ..., 24
C           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
C                    OF DEGREE 12, FOR THE FUNCTION F, IN THE
C                    INTERVAL (BL,BR)
C           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
C                    OF DEGREE 24, FOR THE FUNCTION F, IN THE
C                    INTERVAL (BL,BR)
C           RES12  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB12
C           RES24  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB24
C           QWGTS - EXTERNAL FUNCTION SUBPROGRAM DEFINING
C                    THE FOUR POSSIBLE WEIGHT FUNCTIONS
C           HLGTH  - HALF-LENGTH OF THE INTERVAL (BL,BR)
C           CENTR  - MID POINT OF THE INTERVAL (BL,BR)
C
C***FIRST EXECUTABLE STATEMENT  QC25S
      NEV = 25
      IF(BL.EQ.A.AND.(ALFA.NE.0.0E+00.OR.INTEGR.EQ.2.OR.INTEGR.EQ.4))
     1 GO TO 10
      IF(BR.EQ.B.AND.(BETA.NE.0.0E+00.OR.INTEGR.EQ.3.OR.INTEGR.EQ.4))
     1 GO TO 140
C
C           IF A.GT.BL AND B.LT.BR, APPLY THE 15-POINT GAUSS-KRONROD
C           SCHEME.
C
C
      CALL QK15W(F,QWGTS,A,B,ALFA,BETA,INTEGR,BL,BR,
     1    RESULT,ABSERR,RESABS,RESASC)
      NEV = 15
      GO TO 270
C
C           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF A = BL.
C           ----------------------------------------------------
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F1 = (0.5*(B+B-BR-A)-0.5*(BR-A)*X)**BETA
C                  *F(0.5*(BR-A)*X+0.5*(BR+A))
C
   10 HLGTH = 0.5E+00*(BR-BL)
      CENTR = 0.5E+00*(BR+BL)
      FIX = B-CENTR
      FVAL(1) = 0.5E+00*F(HLGTH+CENTR)*(FIX-HLGTH)**BETA
      FVAL(13) = F(CENTR)*(FIX**BETA)
      FVAL(25) = 0.5E+00*F(CENTR-HLGTH)*(FIX+HLGTH)**BETA
      DO 20 I=2,12
        U = HLGTH*X(I-1)
        ISYM = 26-I
        FVAL(I) = F(U+CENTR)*(FIX-U)**BETA
        FVAL(ISYM) = F(CENTR-U)*(FIX+U)**BETA
   20 CONTINUE
      FACTOR = HLGTH**(ALFA+0.1E+01)
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      RES12 = 0.0E+00
      RES24 = 0.0E+00
      IF(INTEGR.GT.2) GO TO 70
      CALL QCHEB(X,FVAL,CHEB12,CHEB24)
C
C           INTEGR = 1  (OR 2)
C
      DO 30 I=1,13
        RES12 = RES12+CHEB12(I)*RI(I)
        RES24 = RES24+CHEB24(I)*RI(I)
   30 CONTINUE
      DO 40 I=14,25
        RES24 = RES24+CHEB24(I)*RI(I)
   40 CONTINUE
      IF(INTEGR.EQ.1) GO TO 130
C
C           INTEGR = 2
C
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0E+00
      RES24 = 0.0E+00
      DO 50 I=1,13
        RES12 = RES12+CHEB12(I)*RG(I)
        RES24 = RES12+CHEB24(I)*RG(I)
   50 CONTINUE
      DO 60 I=14,25
        RES24 = RES24+CHEB24(I)*RG(I)
   60 CONTINUE
      GO TO 130
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F4 = F1*LOG(0.5*(B+B-BR-A)-0.5*(BR-A)*X)
C
   70 FVAL(1) = FVAL(1)*LOG(FIX-HLGTH)
      FVAL(13) = FVAL(13)*LOG(FIX)
      FVAL(25) = FVAL(25)*LOG(FIX+HLGTH)
      DO 80 I=2,12
        U = HLGTH*X(I-1)
        ISYM = 26-I
        FVAL(I) = FVAL(I)*LOG(FIX-U)
        FVAL(ISYM) = FVAL(ISYM)*LOG(FIX+U)
   80 CONTINUE
      CALL QCHEB(X,FVAL,CHEB12,CHEB24)
C
C           INTEGR = 3  (OR 4)
C
      DO 90 I=1,13
        RES12 = RES12+CHEB12(I)*RI(I)
        RES24 = RES24+CHEB24(I)*RI(I)
   90 CONTINUE
      DO 100 I=14,25
        RES24 = RES24+CHEB24(I)*RI(I)
  100 CONTINUE
      IF(INTEGR.EQ.3) GO TO 130
C
C           INTEGR = 4
C
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0E+00
      RES24 = 0.0E+00
      DO 110 I=1,13
        RES12 = RES12+CHEB12(I)*RG(I)
        RES24 = RES24+CHEB24(I)*RG(I)
  110 CONTINUE
      DO 120 I=14,25
        RES24 = RES24+CHEB24(I)*RG(I)
  120 CONTINUE
  130 RESULT = (RESULT+RES24)*FACTOR
      ABSERR = (ABSERR+ABS(RES24-RES12))*FACTOR
      GO TO 270
C
C           THIS PART OF THE PROGRAM IS EXECUTED ONLY IF B = BR.
C           ----------------------------------------------------
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F2 = (0.5*(B+BL-A-A)+0.5*(B-BL)*X)**ALFA
C                *F(0.5*(B-BL)*X+0.5*(B+BL))
C
  140 HLGTH = 0.5E+00*(BR-BL)
      CENTR = 0.5E+00*(BR+BL)
      FIX = CENTR-A
      FVAL(1) = 0.5E+00*F(HLGTH+CENTR)*(FIX+HLGTH)**ALFA
      FVAL(13) = F(CENTR)*(FIX**ALFA)
      FVAL(25) = 0.5E+00*F(CENTR-HLGTH)*(FIX-HLGTH)**ALFA
      DO 150 I=2,12
        U = HLGTH*X(I-1)
        ISYM = 26-I
        FVAL(I) = F(U+CENTR)*(FIX+U)**ALFA
        FVAL(ISYM) = F(CENTR-U)*(FIX-U)**ALFA
  150 CONTINUE
      FACTOR = HLGTH**(BETA+0.1E+01)
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      RES12 = 0.0E+00
      RES24 = 0.0E+00
      IF(INTEGR.EQ.2.OR.INTEGR.EQ.4) GO TO 200
C
C           INTEGR = 1  (OR 3)
C
      CALL QCHEB(X,FVAL,CHEB12,CHEB24)
      DO 160 I=1,13
        RES12 = RES12+CHEB12(I)*RJ(I)
        RES24 = RES24+CHEB24(I)*RJ(I)
  160 CONTINUE
      DO 170 I=14,25
        RES24 = RES24+CHEB24(I)*RJ(I)
  170 CONTINUE
      IF(INTEGR.EQ.1) GO TO 260
C
C           INTEGR = 3
C
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0E+00
      RES24 = 0.0E+00
      DO 180 I=1,13
        RES12 = RES12+CHEB12(I)*RH(I)
        RES24 = RES24+CHEB24(I)*RH(I)
  180 CONTINUE
      DO 190 I=14,25
        RES24 = RES24+CHEB24(I)*RH(I)
  190 CONTINUE
      GO TO 260
C
C           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
C           FOLLOWING FUNCTION
C           F3 = F2*LOG(0.5*(B-BL)*X+0.5*(B+BL-A-A))
C
  200 FVAL(1) = FVAL(1)*LOG(FIX+HLGTH)
      FVAL(13) = FVAL(13)*LOG(FIX)
      FVAL(25) = FVAL(25)*LOG(FIX-HLGTH)
      DO 210 I=2,12
        U = HLGTH*X(I-1)
        ISYM = 26-I
        FVAL(I) = FVAL(I)*LOG(FIX+U)
        FVAL(ISYM) = FVAL(ISYM)*LOG(FIX-U)
  210 CONTINUE
      CALL QCHEB(X,FVAL,CHEB12,CHEB24)
C
C           INTEGR = 2  (OR 4)
C
      DO 220 I=1,13
        RES12 = RES12+CHEB12(I)*RJ(I)
        RES24 = RES24+CHEB24(I)*RJ(I)
  220 CONTINUE
      DO 230 I=14,25
        RES24 = RES24+CHEB24(I)*RJ(I)
  230 CONTINUE
      IF(INTEGR.EQ.2) GO TO 260
      DC = LOG(BR-BL)
      RESULT = RES24*DC
      ABSERR = ABS((RES24-RES12)*DC)
      RES12 = 0.0E+00
      RES24 = 0.0E+00
C
C           INTEGR = 4
C
      DO 240 I=1,13
        RES12 = RES12+CHEB12(I)*RH(I)
        RES24 = RES24+CHEB24(I)*RH(I)
  240 CONTINUE
      DO 250 I=14,25
        RES24 = RES24+CHEB24(I)*RH(I)
  250 CONTINUE
  260 RESULT = (RESULT+RES24)*FACTOR
      ABSERR = (ABSERR+ABS(RES24-RES12))*FACTOR
  270 RETURN
      END
