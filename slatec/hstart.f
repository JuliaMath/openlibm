*DECK HSTART
      SUBROUTINE HSTART (F, NEQ, A, B, Y, YPRIME, ETOL, MORDER, SMALL,
     +   BIG, SPY, PV, YP, SF, RPAR, IPAR, H)
C***BEGIN PROLOGUE  HSTART
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DEABM, DEBDF and DERKF
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (HSTART-S, DHSTRT-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C   HSTART computes a starting step size to be used in solving initial
C   value problems in ordinary differential equations.
C **********************************************************************
C  Abstract
C
C     Subroutine HSTART computes a starting step size to be used by an
C     initial value method in solving ordinary differential equations.
C     It is based on an estimate of the local Lipschitz constant for the
C     differential equation (lower bound on a norm of the Jacobian),
C     a bound on the differential equation (first derivative), and
C     a bound on the partial derivative of the equation with respect to
C     the independent variable.
C     (All approximated near the initial point A.)
C
C     Subroutine HSTART uses a function subprogram HVNRM for computing
C     a vector norm.  The maximum norm is presently utilized though it
C     can easily be replaced by any other vector norm.  It is presumed
C     that any replacement norm routine would be carefully coded to
C     prevent unnecessary underflows or overflows from occurring, and
C     also, would not alter the vector or number of components.
C
C **********************************************************************
C  On Input you must provide the following
C
C      F -- This is a subroutine of the form
C                               F(X,U,UPRIME,RPAR,IPAR)
C             which defines the system of first order differential
C             equations to be solved.  For the given values of X and the
C             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
C             evaluate the NEQ components of the system of differential
C             equations  dU/DX=F(X,U)  and store the derivatives in the
C             array UPRIME(*), that is,  UPRIME(I) = * dU(I)/DX *  for
C             equations I=1,...,NEQ.
C
C             Subroutine F must not alter X or U(*).  You must declare
C             the name F in an EXTERNAL statement in your program that
C             calls HSTART.  You must dimension U and UPRIME in F.
C
C             RPAR and IPAR are real and integer parameter arrays which
C             you can use for communication between your program and
C             subroutine F.  They are not used or altered by HSTART.  If
C             you do not need RPAR or IPAR, ignore these parameters by
C             treating them as dummy arguments.  If you do choose to use
C             them, dimension them in your program and in F as arrays
C             of appropriate length.
C
C      NEQ -- This is the number of (first order) differential equations
C             to be integrated.
C
C      A -- This is the initial point of integration.
C
C      B -- This is a value of the independent variable used to define
C             the direction of integration.  A reasonable choice is to
C             set  B  to the first point at which a solution is desired.
C             You can also use  B, if necessary, to restrict the length
C             of the first integration step because the algorithm will
C             not compute a starting step length which is bigger than
C             ABS(B-A), unless  B  has been chosen too close to  A.
C             (It is presumed that HSTART has been called with  B
C             different from  A  on the machine being used.  Also see
C             the discussion about the parameter  SMALL.)
C
C      Y(*) -- This is the vector of initial values of the NEQ solution
C             components at the initial point  A.
C
C      YPRIME(*) -- This is the vector of derivatives of the NEQ
C             solution components at the initial point  A.
C             (defined by the differential equations in subroutine F)
C
C      ETOL -- This is the vector of error tolerances corresponding to
C             the NEQ solution components.  It is assumed that all
C             elements are positive.  Following the first integration
C             step, the tolerances are expected to be used by the
C             integrator in an error test which roughly requires that
C                        ABS(local error) .LE. ETOL
C             for each vector component.
C
C      MORDER -- This is the order of the formula which will be used by
C             the initial value method for taking the first integration
C             step.
C
C      SMALL -- This is a small positive machine dependent constant
C             which is used for protecting against computations with
C             numbers which are too small relative to the precision of
C             floating point arithmetic.  SMALL  should be set to
C             (approximately) the smallest positive real number such
C             that  (1.+SMALL) .GT. 1.  on the machine being used. the
C             quantity  SMALL**(3/8)  is used in computing increments of
C             variables for approximating derivatives by differences.
C             also the algorithm will not compute a starting step length
C             which is smaller than  100*SMALL*ABS(A).
C
C      BIG -- This is a large positive machine dependent constant which
C             is used for preventing machine overflows.  A reasonable
C             choice is to set big to (approximately) the square root of
C             the largest real number which can be held in the machine.
C
C      SPY(*),PV(*),YP(*),SF(*) -- These are real work arrays of length
C             NEQ which provide the routine with needed storage space.
C
C      RPAR,IPAR -- These are parameter arrays, of real and integer
C             type, respectively, which can be used for communication
C             between your program and the F subroutine.  They are not
C             used or altered by HSTART.
C
C **********************************************************************
C  On Output  (after the return from HSTART),
C
C      H -- Is an appropriate starting step size to be attempted by the
C             differential equation method.
C
C           All parameters in the call list remain unchanged except for
C           the working arrays SPY(*),PV(*),YP(*) and SF(*).
C
C **********************************************************************
C
C***SEE ALSO  DEABM, DEBDF, DERKF
C***ROUTINES CALLED  HVNRM
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891024  Changed references from VNORM to HVNRM.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  HSTART
C
      DIMENSION Y(*),YPRIME(*),ETOL(*),SPY(*),PV(*),YP(*),SF(*),
     1   RPAR(*),IPAR(*)
      EXTERNAL F
C
C.......................................................................
C
C***FIRST EXECUTABLE STATEMENT  HSTART
      DX = B - A
      ABSDX = ABS(DX)
      RELPER = SMALL**0.375
      YNORM = HVNRM(Y,NEQ)
C
C.......................................................................
C
C     COMPUTE A WEIGHTED APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
C     DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
C     INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW. ALSO
C     COMPUTE A WEIGHTED BOUND (FBND) ON THE FIRST DERIVATIVE LOCALLY.
C
      DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),100.*SMALL*ABS(A)),DX)
      IF (DA .EQ. 0.) DA = RELPER*DX
      CALL F(A+DA,Y,SF,RPAR,IPAR)
C
      IF (MORDER .EQ. 1) GO TO 20
      POWER = 2./(MORDER+1)
      DO 10 J=1,NEQ
        WTJ = ETOL(J)**POWER
        SPY(J) = SF(J)/WTJ
        YP(J) = YPRIME(J)/WTJ
   10   PV(J) = SPY(J) - YP(J)
      GO TO 40
C
   20 DO 30 J=1,NEQ
        SPY(J) = SF(J)/ETOL(J)
        YP(J) = YPRIME(J)/ETOL(J)
   30   PV(J) = SPY(J) - YP(J)
C
   40 DELF = HVNRM(PV,NEQ)
      DFDXB = BIG
      IF (DELF .LT. BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
      YPNORM = HVNRM(YP,NEQ)
      FBND = MAX(HVNRM(SPY,NEQ),YPNORM)
C
C.......................................................................
C
C     COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ CONSTANT FOR
C     THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS ALSO REPRESENTS AN
C     ESTIMATE OF THE NORM OF THE JACOBIAN LOCALLY.
C     THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO ESTIMATE THE
C     LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES. THE FIRST
C     PERTURBATION VECTOR IS BASED ON THE INITIAL DERIVATIVES AND
C     DIRECTION OF INTEGRATION. THE SECOND PERTURBATION VECTOR IS
C     FORMED USING ANOTHER EVALUATION OF THE DIFFERENTIAL EQUATION.
C     THE THIRD PERTURBATION VECTOR IS FORMED USING PERTURBATIONS BASED
C     ONLY ON THE INITIAL VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS
C     CHANGED TO NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
C     INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT COMPONENTS
C     OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE CONSISTENT WITH
C     THE SLOPES OF LOCAL SOLUTION CURVES.
C     ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST DERIVATIVE.
C     NO ATTEMPT IS MADE TO KEEP THE PERTURBATION VECTOR SIZE CONSTANT.
C
      IF (YPNORM .EQ. 0.) GO TO 60
C                       USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
      ICASE = 1
      DO 50 J=1,NEQ
        SPY(J) = YPRIME(J)
   50   YP(J) = YPRIME(J)
      GO TO 80
C                       CANNOT HAVE A NULL PERTURBATION VECTOR
   60 ICASE = 2
      DO 70 J=1,NEQ
        SPY(J) = YPRIME(J)
   70   YP(J) = ETOL(J)
C
   80 DFDUB = 0.
      LK = MIN(NEQ+1,3)
      DO 260 K=1,LK
C                       SET YPNORM AND DELX
        YPNORM = HVNRM(YP,NEQ)
        IF (ICASE .EQ. 1  .OR.  ICASE .EQ. 3) GO TO 90
        DELX = SIGN(1.0,DX)
        GO TO 120
C                       TRY TO ENFORCE MEANINGFUL PERTURBATION VALUES
   90   DELX = DX
        IF (ABS(DELX)*YPNORM .GE. RELPER*YNORM) GO TO 100
        DELXB = BIG
        IF (RELPER*YNORM .LT. BIG*YPNORM) DELXB = RELPER*YNORM/YPNORM
        DELX = SIGN(DELXB,DX)
  100   DO 110 J=1,NEQ
          IF (ABS(DELX*YP(J)) .GT. ETOL(J)) DELX=SIGN(ETOL(J)/YP(J),DX)
  110     CONTINUE
C                       DEFINE PERTURBED VECTOR OF INITIAL VALUES
  120   DO 130 J=1,NEQ
  130     PV(J) = Y(J) + DELX*YP(J)
        IF (K .EQ. 2) GO TO 150
C                       EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
C                       VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
        CALL F(A,PV,YP,RPAR,IPAR)
        DO 140 J=1,NEQ
  140     PV(J) = YP(J) - YPRIME(J)
        GO TO 170
C                       USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
C                                             IN COMPUTING ONE ESTIMATE
  150   CALL F(A+DA,PV,YP,RPAR,IPAR)
        DO 160 J=1,NEQ
  160     PV(J) = YP(J) - SF(J)
C                       CHOOSE LARGEST BOUND ON THE WEIGHTED FIRST
C                                                   DERIVATIVE
  170   IF (MORDER .EQ. 1) GO TO 190
        DO 180 J=1,NEQ
  180     YP(J) = YP(J)/ETOL(J)**POWER
        GO TO 210
  190   DO 200 J=1,NEQ
  200     YP(J) = YP(J)/ETOL(J)
  210   FBND = MAX(FBND,HVNRM(YP,NEQ))
C                       COMPUTE BOUND ON A LOCAL LIPSCHITZ CONSTANT
        DELF = HVNRM(PV,NEQ)
        IF (DELF .EQ. 0.) GO TO 220
        DELY = ABS(DELX)*YPNORM
        IF (DELF .GE. BIG*DELY) GO TO 270
        DFDUB = MAX(DFDUB,DELF/DELY)
C
  220   IF (K .EQ. LK) GO TO 280
C                       CHOOSE NEXT PERTURBATION VECTOR
        DO 250 J=1,NEQ
          IF (K .EQ. LK-1) GO TO 230
          ICASE = 3
          DY = ABS(PV(J))
          IF (DY .EQ. 0.) DY = MAX(DELF,ETOL(J))
          GO TO 240
  230     ICASE = 4
          DY = MAX(RELPER*ABS(Y(J)),ETOL(J))
  240     IF (SPY(J) .EQ. 0.) SPY(J) = YP(J)
          IF (SPY(J) .NE. 0.) DY = SIGN(DY,SPY(J))
  250     YP(J) = DY
  260   CONTINUE
C
C                       PROTECT AGAINST AN OVERFLOW
  270 DFDUB = BIG
C
C.......................................................................
C
C     COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
C
  280 YDPB = DFDXB + DFDUB*FBND
C
C.......................................................................
C
C     COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND SECOND
C     DERIVATIVE INFORMATION
C
C                       RESTRICT THE STEP LENGTH TO BE NOT BIGGER THAN
C                       ABS(B-A).   (UNLESS  B  IS TOO CLOSE TO  A)
      H = ABSDX
C
      IF (YDPB .NE. 0.  .OR.  FBND .NE. 0.) GO TO 290
C
C                       BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
C                                    DERIVATIVE TERM (YDPB) ARE ZERO
      GO TO 310
C
  290 IF (YDPB .NE. 0.) GO TO 300
C
C                       ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
      IF (1.0 .LT. FBND*ABSDX) H = 1./FBND
      GO TO 310
C
C                       SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
  300 SRYDPB = SQRT(0.5*YDPB)
      IF (1.0 .LT. SRYDPB*ABSDX) H = 1./SRYDPB
C
C                       FURTHER RESTRICT THE STEP LENGTH TO BE NOT
C                                                 BIGGER THAN  1/DFDUB
  310 IF (H*DFDUB .GT. 1.) H = 1./DFDUB
C
C                       FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
C                       SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
C                       A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
C                       THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
C                                                       STEP LENGTH.
      H = MAX(H,100.*SMALL*ABS(A))
      IF (H .EQ. 0.) H = SMALL*ABS(B)
C
C                       NOW SET DIRECTION OF INTEGRATION
      H = SIGN(H,DX)
C
      RETURN
      END
