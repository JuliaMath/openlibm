*DECK STEPS
      SUBROUTINE STEPS (F, NEQN, Y, X, H, EPS, WT, START, HOLD, K, KOLD,
     +   CRASH, PHI, P, YP, PSI, ALPHA, BETA, SIG, V, W, G, PHASE1, NS,
     +   NORND, KSTEPS, TWOU, FOURU, XOLD, KPREV, IVC, IV, KGI, GI,
     +   RPAR, IPAR)
C***BEGIN PROLOGUE  STEPS
C***PURPOSE  Integrate a system of first order ordinary differential
C            equations one step.
C***LIBRARY   SLATEC (DEPAC)
C***CATEGORY  I1A1B
C***TYPE      SINGLE PRECISION (STEPS-S, DSTEPS-D)
C***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
C             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
C***AUTHOR  Shampine, L. F., (SNLA)
C           Gordon, M. K., (SNLA)
C             MODIFIED BY H.A. WATTS
C***DESCRIPTION
C
C   Written by L. F. Shampine and M. K. Gordon
C
C   Abstract
C
C   Subroutine  STEPS  is normally used indirectly through subroutine
C   DEABM .  Because  DEABM  suffices for most problems and is much
C   easier to use, using it should be considered before using  STEPS
C   alone.
C
C   Subroutine STEPS integrates a system of  NEQN  first order ordinary
C   differential equations one step, normally from X to X+H, using a
C   modified divided difference form of the Adams Pece formulas.  Local
C   extrapolation is used to improve absolute stability and accuracy.
C   The code adjusts its order and step size to control the local error
C   per unit step in a generalized sense.  Special devices are included
C   to control roundoff error and to detect when the user is requesting
C   too much accuracy.
C
C   This code is completely explained and documented in the text,
C   Computer Solution of Ordinary Differential Equations, The Initial
C   Value Problem  by L. F. Shampine and M. K. Gordon.
C   Further details on use of this code are available in "Solving
C   Ordinary Differential Equations with ODE, STEP, and INTRP",
C   by L. F. Shampine and M. K. Gordon, SLA-73-1060.
C
C
C   The parameters represent --
C      F -- subroutine to evaluate derivatives
C      NEQN -- number of equations to be integrated
C      Y(*) -- solution vector at X
C      X -- independent variable
C      H -- appropriate step size for next step.  Normally determined by
C           code
C      EPS -- local error tolerance
C      WT(*) -- vector of weights for error criterion
C      START -- logical variable set .TRUE. for first step,  .FALSE.
C           otherwise
C      HOLD -- step size used for last successful step
C      K -- appropriate order for next step (determined by code)
C      KOLD -- order used for last successful step
C      CRASH -- logical variable set .TRUE. when no step can be taken,
C           .FALSE. otherwise.
C      YP(*) -- derivative of solution vector at  X  after successful
C           step
C      KSTEPS -- counter on attempted steps
C      TWOU -- 2.*U where U is machine unit roundoff quantity
C      FOURU -- 4.*U where U is machine unit roundoff quantity
C      RPAR,IPAR -- parameter arrays which you may choose to use
C            for communication between your program and subroutine F.
C            They are not altered or used by STEPS.
C   The variables X,XOLD,KOLD,KGI and IVC and the arrays Y,PHI,ALPHA,G,
C   W,P,IV and GI are required for the interpolation subroutine SINTRP.
C   The remaining variables and arrays are included in the call list
C   only to eliminate local retention of variables between calls.
C
C   Input to STEPS
C
C      First call --
C
C   The user must provide storage in his calling program for all arrays
C   in the call list, namely
C
C     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12),
C    1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
C    2  RPAR(*),IPAR(*)
C
C    **Note**
C
C   The user must also declare  START ,  CRASH ,  PHASE1  and  NORND
C   logical variables and  F  an EXTERNAL subroutine, supply the
C   subroutine  F(X,Y,YP)  to evaluate
C      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN))
C   and initialize only the following parameters.
C      NEQN -- number of equations to be integrated
C      Y(*) -- vector of initial values of dependent variables
C      X -- initial value of the independent variable
C      H -- nominal step size indicating direction of integration
C           and maximum size of step.  Must be variable
C      EPS -- local error tolerance per step.  Must be variable
C      WT(*) -- vector of non-zero weights for error criterion
C      START -- .TRUE.
C      YP(*) -- vector of initial derivative values
C      KSTEPS -- set KSTEPS to zero
C      TWOU -- 2.*U where U is machine unit roundoff quantity
C      FOURU -- 4.*U where U is machine unit roundoff quantity
C   Define U to be the machine unit roundoff quantity by calling
C   the function routine  R1MACH,  U = R1MACH(4), or by
C   computing U so that U is the smallest positive number such
C   that 1.0+U .GT. 1.0.
C
C   STEPS  requires that the L2 norm of the vector with components
C   LOCAL ERROR(L)/WT(L)  be less than  EPS  for a successful step.  The
C   array  WT  allows the user to specify an error test appropriate
C   for his problem.  For example,
C      WT(L) = 1.0  specifies absolute error,
C            = ABS(Y(L))  error relative to the most recent value of the
C                 L-th component of the solution,
C            = ABS(YP(L))  error relative to the most recent value of
C                 the L-th component of the derivative,
C            = MAX(WT(L),ABS(Y(L)))  error relative to the largest
C                 magnitude of L-th component obtained so far,
C            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  specifies a mixed
C                 relative-absolute test where  RELERR  is relative
C                 error,  ABSERR  is absolute error and  EPS =
C                 MAX(RELERR,ABSERR) .
C
C      Subsequent calls --
C
C   Subroutine  STEPS  is designed so that all information needed to
C   continue the integration, including the step size  H  and the order
C   K , is returned with each step.  With the exception of the step
C   size, the error tolerance, and the weights, none of the parameters
C   should be altered.  The array  WT  must be updated after each step
C   to maintain relative error tests like those above.  Normally the
C   integration is continued just beyond the desired endpoint and the
C   solution interpolated there with subroutine  SINTRP .  If it is
C   impossible to integrate beyond the endpoint, the step size may be
C   reduced to hit the endpoint since the code will not take a step
C   larger than the  H  input.  Changing the direction of integration,
C   i.e., the sign of  H , requires the user set  START = .TRUE. before
C   calling  STEPS  again.  This is the only situation in which  START
C   should be altered.
C
C   Output from STEPS
C
C      Successful Step --
C
C   The subroutine returns after each successful step with  START  and
C   CRASH  set .FALSE. .  X  represents the independent variable
C   advanced one step of length  HOLD  from its value on input and  Y
C   the solution vector at the new value of  X .  All other parameters
C   represent information corresponding to the new  X  needed to
C   continue the integration.
C
C      Unsuccessful Step --
C
C   When the error tolerance is too small for the machine precision,
C   the subroutine returns without taking a step and  CRASH = .TRUE. .
C   An appropriate step size and error tolerance for continuing are
C   estimated and all other information is restored as upon input
C   before returning.  To continue with the larger tolerance, the user
C   just calls the code again.  A restart is neither required nor
C   desirable.
C
C***REFERENCES  L. F. Shampine and M. K. Gordon, Solving ordinary
C                 differential equations with ODE, STEP, and INTRP,
C                 Report SLA-73-1060, Sandia Laboratories, 1973.
C***ROUTINES CALLED  HSTART, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   740101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  STEPS
C
      LOGICAL START,CRASH,PHASE1,NORND
      DIMENSION Y(*),WT(*),PHI(NEQN,16),P(*),YP(*),PSI(12),
     1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
     2  RPAR(*),IPAR(*)
      DIMENSION TWO(13),GSTR(13)
      EXTERNAL F
      SAVE TWO, GSTR
C
      DATA  TWO(1),TWO(2),TWO(3),TWO(4),TWO(5),TWO(6),TWO(7),TWO(8),
     1      TWO(9),TWO(10),TWO(11),TWO(12),TWO(13) /2.0,4.0,8.0,16.0,
     2      32.0,64.0,128.0,256.0,512.0,1024.0,2048.0,4096.0,8192.0/
      DATA  GSTR(1),GSTR(2),GSTR(3),GSTR(4),GSTR(5),GSTR(6),GSTR(7),
     1      GSTR(8),GSTR(9),GSTR(10),GSTR(11),GSTR(12),GSTR(13)/0.500,
     2      0.0833,0.0417,0.0264,0.0188,0.0143,0.0114,0.00936,0.00789,
     3      0.00679,0.00592,0.00524,0.00468/
C
C
C       ***     BEGIN BLOCK 0     ***
C   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
C   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
C   STARTING STEP SIZE.
C                   ***
C
C   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
C
C***FIRST EXECUTABLE STATEMENT  STEPS
      CRASH = .TRUE.
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5*EPS
C
C   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE
C
      ROUND = 0.0
      DO 10 L = 1,NEQN
 10     ROUND = ROUND + (Y(L)/WT(L))**2
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS .GE. ROUND) GO TO 15
      EPS = 2.0*ROUND*(1.0 + FOURU)
      RETURN
 15   CRASH = .FALSE.
      G(1) = 1.0
      G(2) = 0.5
      SIG(1) = 1.0
      IF(.NOT.START) GO TO 99
C
C   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP
C
C     CALL F(X,Y,YP,RPAR,IPAR)
C     SUM = 0.0
      DO 20 L = 1,NEQN
        PHI(L,1) = YP(L)
   20   PHI(L,2) = 0.0
C20     SUM = SUM + (YP(L)/WT(L))**2
C     SUM = SQRT(SUM)
C     ABSH = ABS(H)
C     IF(EPS .LT. 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM)
C     H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
C
      U = R1MACH(4)
      BIG = SQRT(R1MACH(2))
      CALL HSTART (F,NEQN,X,X+H,Y,YP,WT,1,U,BIG,
     1             PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H)
C
      HOLD = 0.0
      K = 1
      KOLD = 0
      KPREV = 0
      START = .FALSE.
      PHASE1 = .TRUE.
      NORND = .TRUE.
      IF(P5EPS .GT. 100.0*ROUND) GO TO 99
      NORND = .FALSE.
      DO 25 L = 1,NEQN
 25     PHI(L,15) = 0.0
 99   IFAIL = 0
C       ***     END BLOCK 0     ***
C
C       ***     BEGIN BLOCK 1     ***
C   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
C   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED.
C                   ***
C
 100  KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      KM2 = K-2
C
C   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
C   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
C
      IF(H .NE. HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1
      NSP1 = NS+1
      IF (K .LT. NS) GO TO 199
C
C   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
C   ARE CHANGED
C
      BETA(NS) = 1.0
      REALNS = NS
      ALPHA(NS) = 1.0/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = 1.0
      IF(K .LT. NSP1) GO TO 110
      DO 105 I = NSP1,K
        IM1 = I-1
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H
        ALPHA(I) = H/TEMP1
        REALI = I
 105    SIG(I+1) = REALI*ALPHA(I)*SIG(I)
 110  PSI(K) = TEMP1
C
C   COMPUTE COEFFICIENTS G(*)
C
C   INITIALIZE V(*) AND SET W(*).
C
      IF(NS .GT. 1) GO TO 120
      DO 115 IQ = 1,K
        TEMP3 = IQ*(IQ+1)
        V(IQ) = 1.0/TEMP3
 115    W(IQ) = V(IQ)
      IVC = 0
      KGI = 0
      IF (K .EQ. 1) GO TO 140
      KGI = 1
      GI(1) = W(2)
      GO TO 140
C
C   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*)
C
 120  IF(K .LE. KPREV) GO TO 130
      IF (IVC .EQ. 0) GO TO 122
      JV = KP1 - IV(IVC)
      IVC = IVC - 1
      GO TO 123
 122  JV = 1
      TEMP4 = K*KP1
      V(K) = 1.0/TEMP4
      W(K) = V(K)
      IF (K .NE. 2) GO TO 123
      KGI = 1
      GI(1) = W(2)
 123  NSM2 = NS-2
      IF(NSM2 .LT. JV) GO TO 130
      DO 125 J = JV,NSM2
        I = K-J
        V(I) = V(I) - ALPHA(J+1)*V(I+1)
 125    W(I) = V(I)
      IF (I .NE. 2) GO TO 130
      KGI = NS - 1
      GI(KGI) = W(2)
C
C   UPDATE V(*) AND SET W(*)
C
 130  LIMIT1 = KP1 - NS
      TEMP5 = ALPHA(NS)
      DO 135 IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1)
 135    W(IQ) = V(IQ)
      G(NSP1) = W(1)
      IF (LIMIT1 .EQ. 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1)
      IF (K .GE. KOLD) GO TO 140
      IVC = IVC + 1
      IV(IVC) = LIMIT1 + 2
C
C   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
C
 140  NSP2 = NS + 2
      KPREV = K
      IF(KP1 .LT. NSP2) GO TO 199
      DO 150 I = NSP2,KP1
        LIMIT2 = KP2 - I
        TEMP6 = ALPHA(I-1)
        DO 145 IQ = 1,LIMIT2
 145      W(IQ) = W(IQ) - TEMP6*W(IQ+1)
 150    G(I) = W(1)
 199    CONTINUE
C       ***     END BLOCK 1     ***
C
C       ***     BEGIN BLOCK 2     ***
C   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED
C   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K,
C   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
C                   ***
C
C   INCREMENT COUNTER ON ATTEMPTED STEPS
C
      KSTEPS = KSTEPS + 1
C
C   CHANGE PHI TO PHI STAR
C
      IF(K .LT. NSP1) GO TO 215
      DO 210 I = NSP1,K
        TEMP1 = BETA(I)
        DO 205 L = 1,NEQN
 205      PHI(L,I) = TEMP1*PHI(L,I)
 210    CONTINUE
C
C   PREDICT SOLUTION AND DIFFERENCES
C
 215  DO 220 L = 1,NEQN
        PHI(L,KP2) = PHI(L,KP1)
        PHI(L,KP1) = 0.0
 220    P(L) = 0.0
      DO 230 J = 1,K
        I = KP1 - J
        IP1 = I+1
        TEMP2 = G(I)
        DO 225 L = 1,NEQN
          P(L) = P(L) + TEMP2*PHI(L,I)
 225      PHI(L,I) = PHI(L,I) + PHI(L,IP1)
 230    CONTINUE
      IF(NORND) GO TO 240
      DO 235 L = 1,NEQN
        TAU = H*P(L) - PHI(L,15)
        P(L) = Y(L) + TAU
 235    PHI(L,16) = (P(L) - Y(L)) - TAU
      GO TO 250
 240  DO 245 L = 1,NEQN
 245    P(L) = Y(L) + H*P(L)
 250  XOLD = X
      X = X + H
      ABSH = ABS(H)
      CALL F(X,P,YP,RPAR,IPAR)
C
C   ESTIMATE ERRORS AT ORDERS K,K-1,K-2
C
      ERKM2 = 0.0
      ERKM1 = 0.0
      ERK = 0.0
      DO 265 L = 1,NEQN
        TEMP3 = 1.0/WT(L)
        TEMP4 = YP(L) - PHI(L,1)
        IF(KM2)265,260,255
 255    ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
 260    ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2
 265    ERK = ERK + (TEMP4*TEMP3)**2
      IF(KM2)280,275,270
 270  ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2)
 275  ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1)
 280  TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1))
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
C
C   TEST IF ORDER SHOULD BE LOWERED
C
      IF(KM2)299,290,285
 285  IF(MAX(ERKM1,ERKM2) .LE. ERK) KNEW = KM1
      GO TO 299
 290  IF(ERKM1 .LE. 0.5*ERK) KNEW = KM1
C
C   TEST IF STEP SUCCESSFUL
C
 299  IF(ERR .LE. EPS) GO TO 400
C       ***     END BLOCK 2     ***
C
C       ***     BEGIN BLOCK 3     ***
C   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) .
C   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE
C   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
C   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
C   PRECISION.
C                   ***
C
C   RESTORE X, PHI(*,*) AND PSI(*)
C
      PHASE1 = .FALSE.
      X = XOLD
      DO 310 I = 1,K
        TEMP1 = 1.0/BETA(I)
        IP1 = I+1
        DO 305 L = 1,NEQN
 305      PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
 310    CONTINUE
      IF(K .LT. 2) GO TO 320
      DO 315 I = 2,K
 315    PSI(I-1) = PSI(I) - H
C
C   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP
C   SIZE
C
 320  IFAIL = IFAIL + 1
      TEMP2 = 0.5
      IF(IFAIL - 3) 335,330,325
 325  IF(P5EPS .LT. 0.25*ERK) TEMP2 = SQRT(P5EPS/ERK)
 330  KNEW = 1
 335  H = TEMP2*H
      K = KNEW
      NS = 0
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS
      RETURN
 340  GO TO 100
C       ***     END BLOCK 3     ***
C
C       ***     BEGIN BLOCK 4     ***
C   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE
C   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE
C   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP.
C                   ***
 400  KOLD = K
      HOLD = H
C
C   CORRECT AND EVALUATE
C
      TEMP1 = H*G(KP1)
      IF(NORND) GO TO 410
      DO 405 L = 1,NEQN
        TEMP3 = Y(L)
        RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
        Y(L) = P(L) + RHO
        PHI(L,15) = (Y(L) - P(L)) - RHO
 405    P(L) = TEMP3
      GO TO 420
 410  DO 415 L = 1,NEQN
        TEMP3 = Y(L)
        Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
 415    P(L) = TEMP3
 420  CALL F(X,Y,YP,RPAR,IPAR)
C
C   UPDATE DIFFERENCES FOR NEXT STEP
C
      DO 425 L = 1,NEQN
        PHI(L,KP1) = YP(L) - PHI(L,1)
 425    PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
      DO 435 I = 1,K
        DO 430 L = 1,NEQN
 430      PHI(L,I) = PHI(L,I) + PHI(L,KP1)
 435    CONTINUE
C
C   ESTIMATE ERROR AT ORDER K+1 UNLESS:
C     IN FIRST PHASE WHEN ALWAYS RAISE ORDER,
C     ALREADY DECIDED TO LOWER ORDER,
C     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE
C
      ERKP1 = 0.0
      IF(KNEW .EQ. KM1  .OR.  K .EQ. 12) PHASE1 = .FALSE.
      IF(PHASE1) GO TO 450
      IF(KNEW .EQ. KM1) GO TO 455
      IF(KP1 .GT. NS) GO TO 460
      DO 440 L = 1,NEQN
 440    ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
C
C   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER
C   FOR NEXT STEP
C
      IF(K .GT. 1) GO TO 445
      IF(ERKP1 .GE. 0.5*ERK) GO TO 460
      GO TO 450
 445  IF(ERKM1 .LE. MIN(ERK,ERKP1)) GO TO 455
      IF(ERKP1 .GE. ERK  .OR.  K .EQ. 12) GO TO 460
C
C   HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE
C   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
C
C   RAISE ORDER
C
 450  K = KP1
      ERK = ERKP1
      GO TO 460
C
C   LOWER ORDER
C
 455  K = KM1
      ERK = ERKM1
C
C   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
C
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS .GE. ERK*TWO(K+1)) GO TO 465
      HNEW = H
      IF(P5EPS .GE. ERK) GO TO 465
      TEMP2 = K+1
      R = (P5EPS/ERK)**(1.0/TEMP2)
      HNEW = ABSH*MAX(0.5,MIN(0.9,R))
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H)
 465  H = HNEW
      RETURN
C       ***     END BLOCK 4     ***
      END
