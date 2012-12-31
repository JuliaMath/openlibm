*DECK SOSEQS
      SUBROUTINE SOSEQS (FNC, N, S, RTOLX, ATOLX, TOLF, IFLAG, MXIT,
     +   NCJS, NSRRC, NSRI, IPRINT, FMAX, C, NC, B, P, TEMP, X, Y, FAC,
     +   IS)
C***BEGIN PROLOGUE  SOSEQS
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SOS
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SOSEQS-S, DSOSEQ-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     SOSEQS solves a system of N simultaneous nonlinear equations.
C     See the comments in the interfacing routine SOS for a more
C     detailed description of some of the items in the calling list.
C
C ********************************************************************
C
C   -INPUT-
C     FNC -Function subprogram which evaluates the equations
C     N   -Number of equations
C     S   -Solution vector of initial guesses
C     RTOLX-Relative error tolerance on solution components
C     ATOLX-Absolute error tolerance on solution components
C     TOLF-Residual error tolerance
C     MXIT-Maximum number of allowable iterations.
C     NCJS-Maximum number of consecutive iterative steps to perform
C          using the same triangular Jacobian matrix approximation.
C     NSRRC-Number of consecutive iterative steps for which the
C          limiting precision accuracy test must be satisfied
C          before the routine exits with IFLAG=4.
C     NSRI-Number of consecutive iterative steps for which the
C          diverging condition test must be satisfied before
C          the routine exits with IFLAG=7.
C     IPRINT-Internal printing parameter.  You must set IPRINT=-1 if you
C          want the intermediate solution iterates and a residual norm
C          to be printed.
C     C   -Internal work array, dimensioned at least N*(N+1)/2.
C     NC  -Dimension of C array. NC  .GE.  N*(N+1)/2.
C     B   -Internal work array, dimensioned N.
C     P   -Internal work array, dimensioned N.
C     TEMP-Internal work array, dimensioned N.
C     X   -Internal work array, dimensioned N.
C     Y   -Internal work array, dimensioned N.
C     FAC -Internal work array, dimensioned N.
C     IS  -Internal work array, dimensioned N.
C
C   -OUTPUT-
C     S   -Solution vector
C     IFLAG-Status indicator flag
C     MXIT-The actual number of iterations performed
C     FMAX-Residual norm
C     C   -Upper unit triangular matrix which approximates the
C          forward triangularization of the full Jacobian matrix.
C          stored in a vector with dimension at least N*(N+1)/2.
C     B   -Contains the residuals (function values) divided
C          by the corresponding components of the P vector
C     P   -Array used to store the partial derivatives. After
C          each iteration P(K) contains the maximal derivative
C          occurring in the K-th reduced equation.
C     TEMP-Array used to store the previous solution iterate.
C     X   -Solution vector. Contains the values achieved on the
C          last iteration loop upon exit from SOS.
C     Y   -Array containing the solution increments.
C     FAC -Array containing factors used in computing numerical
C          derivatives.
C     IS  -Records the pivotal information (column interchanges)
C
C **********************************************************************
C *** Three machine dependent parameters appear in this subroutine.
C
C *** The smallest positive magnitude, zero, is defined by the function
C *** routine R1MACH(1).
C
C *** URO, The computer unit roundoff value, is defined by R1MACH(3) for
C *** machines that round or R1MACH(4) for machines that truncate.
C *** URO is the smallest positive number such that 1.+URO  .GT.  1.
C
C *** The output tape unit number, LOUN, is defined by the function
C *** I1MACH(2).
C **********************************************************************
C
C***SEE ALSO  SOS
C***ROUTINES CALLED  I1MACH, R1MACH, SOSSOL
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  SOSEQS
C
C
      DIMENSION S(*), C(NC), B(*), IS(*), P(*), TEMP(*), X(*), Y(*),
     1          FAC(*)
C
C***FIRST EXECUTABLE STATEMENT  SOSEQS
      URO = R1MACH(4)
      LOUN = I1MACH(2)
      ZERO = R1MACH(1)
      RE = MAX(RTOLX,URO)
      SRURO = SQRT(URO)
C
      IFLAG = 0
      NP1 = N + 1
      ICR = 0
      IC = 0
      ITRY = NCJS
      YN1 = 0.
      YN2 = 0.
      YN3 = 0.
      YNS = 0.
      MIT = 0
      FN1 = 0.
      FN2 = 0.
      FMXS = 0.
C
C     INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND
C     SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE.
C
      DO 10 K=1,N
        IS(K) = K
        X(K) = S(K)
        TEMP(K) = X(K)
   10 CONTINUE
C
C
C    *****************************************
C    **** BEGIN PRINCIPAL ITERATION LOOP  ****
C    *****************************************
C
      DO 330 M=1,MXIT
C
        DO 20 K=1,N
          FAC(K) = SRURO
   20   CONTINUE
C
   30   KN = 1
        FMAX = 0.
C
C
C    ******** BEGIN SUBITERATION LOOP DEFINING THE LINEARIZATION OF EACH
C    ******** EQUATION WHICH RESULTS IN THE CONSTRUCTION OF AN UPPER
C    ******** TRIANGULAR MATRIX APPROXIMATING THE FORWARD
C    ******** TRIANGULARIZATION OF THE FULL JACOBIAN MATRIX
C
        DO 170 K=1,N
          KM1 = K - 1
C
C     BACK-SOLVE A TRIANGULAR LINEAR SYSTEM OBTAINING
C     IMPROVED SOLUTION VALUES FOR K-1 OF THE VARIABLES
C     FROM THE FIRST K-1 EQUATIONS. THESE VARIABLES ARE THEN
C     ELIMINATED FROM THE K-TH EQUATION.
C
          IF (KM1 .EQ. 0) GO TO 50
          CALL SOSSOL(K, N, KM1, Y, C, B, KN)
          DO 40 J=1,KM1
            JS = IS(J)
            X(JS) = TEMP(JS) + Y(J)
   40     CONTINUE
C
C
C     EVALUATE THE K-TH EQUATION AND THE INTERMEDIATE COMPUTATION
C     FOR THE MAX NORM OF THE RESIDUAL VECTOR.
C
   50     F = FNC(X,K)
          FMAX = MAX(FMAX,ABS(F))
C
C     IF WE WISH TO PERFORM SEVERAL ITERATIONS USING A FIXED
C     FACTORIZATION OF AN APPROXIMATE JACOBIAN,WE NEED ONLY
C     UPDATE THE CONSTANT VECTOR.
C
          IF (ITRY .LT. NCJS) GO TO 160
C
C
          IT = 0
C
C     COMPUTE PARTIAL DERIVATIVES THAT ARE REQUIRED IN THE LINEARIZATION
C     OF THE K-TH REDUCED EQUATION
C
          DO 90 J=K,N
            ITEM = IS(J)
            HX = X(ITEM)
            H = FAC(ITEM)*HX
            IF (ABS(H) .LE. ZERO) H = FAC(ITEM)
            X(ITEM) = HX + H
            IF (KM1 .EQ. 0) GO TO 70
            Y(J) = H
            CALL SOSSOL(K, N, J, Y, C, B, KN)
            DO 60 L=1,KM1
              LS = IS(L)
              X(LS) = TEMP(LS) + Y(L)
   60       CONTINUE
   70       FP = FNC(X,K)
            X(ITEM) = HX
            FDIF = FP - F
            IF (ABS(FDIF) .GT. URO*ABS(F)) GO TO 80
            FDIF = 0.
            IT = IT + 1
   80       P(J) = FDIF/H
   90     CONTINUE
C
          IF (IT .LE. (N-K)) GO TO 110
C
C     ALL COMPUTED PARTIAL DERIVATIVES OF THE K-TH EQUATION
C     ARE EFFECTIVELY ZERO.TRY LARGER PERTURBATIONS OF THE
C     INDEPENDENT VARIABLES.
C
          DO 100 J=K,N
            ISJ = IS(J)
            FACT = 100.*FAC(ISJ)
            IF (FACT .GT. 1.E+10) GO TO 340
            FAC(ISJ) = FACT
  100     CONTINUE
          GO TO 30
C
  110     IF (K .EQ. N) GO TO 160
C
C     ACHIEVE A PIVOTING EFFECT BY CHOOSING THE MAXIMAL DERIVATIVE
C     ELEMENT
C
          PMAX = 0.
          DO 120 J=K,N
            TEST = ABS(P(J))
            IF (TEST .LE. PMAX) GO TO 120
            PMAX = TEST
            ISV = J
  120     CONTINUE
          IF (PMAX .EQ. 0.) GO TO 340
C
C     SET UP THE COEFFICIENTS FOR THE K-TH ROW OF THE TRIANGULAR
C     LINEAR SYSTEM AND SAVE THE PARTIAL DERIVATIVE OF
C     LARGEST MAGNITUDE
C
          PMAX = P(ISV)
          KK = KN
          DO 140 J=K,N
            IF (J .EQ. ISV) GO TO 130
            C(KK) = -P(J)/PMAX
  130       KK = KK + 1
  140     CONTINUE
          P(K) = PMAX
C
C
          IF (ISV .EQ. K) GO TO 160
C
C     INTERCHANGE THE TWO COLUMNS OF C DETERMINED BY THE
C     PIVOTAL STRATEGY
C
          KSV = IS(K)
          IS(K) = IS(ISV)
          IS(ISV) = KSV
C
          KD = ISV - K
          KJ = K
          DO 150 J=1,K
            CSV = C(KJ)
            JK = KJ + KD
            C(KJ) = C(JK)
            C(JK) = CSV
            KJ = KJ + N - J
  150     CONTINUE
C
  160     KN = KN + NP1 - K
C
C     STORE THE COMPONENTS FOR THE CONSTANT VECTOR
C
          B(K) = -F/P(K)
C
  170   CONTINUE
C
C    ********
C    ******** END OF LOOP CREATING THE TRIANGULAR LINEARIZATION MATRIX
C    ********
C
C
C     SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW SOLUTION
C     APPROXIMATION AND OBTAIN THE SOLUTION INCREMENT NORM.
C
        KN = KN - 1
        Y(N) = B(N)
        IF (N .GT. 1) CALL SOSSOL(N, N, N, Y, C, B, KN)
        XNORM = 0.
        YNORM = 0.
        DO 180 J=1,N
          YJ = Y(J)
          YNORM = MAX(YNORM,ABS(YJ))
          JS = IS(J)
          X(JS) = TEMP(JS) + YJ
          XNORM = MAX(XNORM,ABS(X(JS)))
  180   CONTINUE
C
C
C     PRINT INTERMEDIATE SOLUTION ITERATES AND RESIDUAL NORM IF DESIRED
C
        IF (IPRINT.NE.(-1)) GO TO 190
        MM = M - 1
        WRITE (LOUN,1234) FMAX, MM, (X(J),J=1,N)
 1234   FORMAT ('0RESIDUAL NORM =', E9.2, /1X, 'SOLUTION ITERATE',
     1   ' (', I3, ')', /(1X, 5E26.14))
  190   CONTINUE
C
C     TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE AND/OR ABSOLUTE ERROR
C     COMPARISON ON SUCCESSIVE APPROXIMATIONS OF EACH SOLUTION VARIABLE)
C
        DO 200 J=1,N
          JS = IS(J)
          IF (ABS(Y(J)) .GT. RE*ABS(X(JS))+ATOLX) GO TO 210
  200   CONTINUE
        IF (FMAX .LE. FMXS) IFLAG = 1
C
C     TEST FOR CONVERGENCE TO A SOLUTION BASED ON RESIDUALS
C
  210   IF (FMAX .GT. TOLF) GO TO 220
        IFLAG = IFLAG + 2
  220   IF (IFLAG .GT. 0) GO TO 360
C
C
        IF (M .GT. 1) GO TO 230
        FMIN = FMAX
        GO TO 280
C
C     SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM.
C
  230   IF (FMAX .GE. FMIN) GO TO 250
        MIT = M + 1
        YN1 = YNORM
        YN2 = YNS
        FN1 = FMXS
        FMIN = FMAX
        DO 240 J=1,N
          S(J) = X(J)
  240   CONTINUE
        IC = 0
C
C     TEST FOR LIMITING PRECISION CONVERGENCE.  VERY SLOWLY CONVERGENT
C     PROBLEMS MAY ALSO BE DETECTED.
C
  250   IF (YNORM .GT. SRURO*XNORM) GO TO 260
        IF ((FMAX .LT. 0.2*FMXS) .OR. (FMAX .GT. 5.*FMXS)) GO TO 260
        IF ((YNORM .LT. 0.2*YNS) .OR. (YNORM .GT. 5.*YNS)) GO TO 260
        ICR = ICR + 1
        IF (ICR .LT. NSRRC) GO TO 270
        IFLAG = 4
        FMAX = FMIN
        GO TO 380
  260   ICR = 0
C
C     TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME.
C
        IF ((YNORM .LE. 2.*YNS) .AND. (FMAX .LE. 2.*FMXS)) GO TO 270
        IC = IC + 1
        IF (IC .LT. NSRI) GO TO 280
        IFLAG = 7
        GO TO 360
  270   IC = 0
C
C     CHECK TO SEE IF NEXT ITERATION CAN USE THE OLD JACOBIAN
C     FACTORIZATION
C
  280   ITRY = ITRY - 1
        IF (ITRY .EQ. 0) GO TO 290
        IF (20.*YNORM .GT. XNORM) GO TO 290
        IF (YNORM .GT. 2.*YNS) GO TO 290
        IF (FMAX .LT. 2.*FMXS) GO TO 300
  290   ITRY = NCJS
C
C     SAVE THE CURRENT SOLUTION APPROXIMATION AND THE RESIDUAL AND
C     SOLUTION INCREMENT NORMS FOR USE IN THE NEXT ITERATION.
C
  300   DO 310 J=1,N
          TEMP(J) = X(J)
  310   CONTINUE
        IF (M.NE.MIT) GO TO 320
        FN2 = FMAX
        YN3 = YNORM
  320   FMXS = FMAX
        YNS = YNORM
C
C
  330 CONTINUE
C
C    *****************************************
C    **** END OF PRINCIPAL ITERATION LOOP ****
C    *****************************************
C
C
C     TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED.
      M = MXIT
      IFLAG = 5
      IF (YN1 .GT. 10.0*YN2 .OR. YN3 .GT. 10.0*YN1) IFLAG = 6
      IF (FN1 .GT. 5.0*FMIN .OR. FN2 .GT. 5.0*FMIN) IFLAG = 6
      IF (FMAX .GT. 5.0*FMIN) IFLAG = 6
      GO TO 360
C
C
C     A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR.
  340 IFLAG = 8
      DO 350 J=1,N
        S(J) = TEMP(J)
  350 CONTINUE
      GO TO 380
C
C
  360 DO 370 J=1,N
        S(J) = X(J)
  370 CONTINUE
C
C
  380 MXIT = M
      RETURN
      END
