*DECK DSOSEQ
      SUBROUTINE DSOSEQ (FNC, N, S, RTOLX, ATOLX, TOLF, IFLAG, MXIT,
     +   NCJS, NSRRC, NSRI, IPRINT, FMAX, C, NC, B, P, TEMP, X, Y, FAC,
     +   IS)
C***BEGIN PROLOGUE  DSOSEQ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DSOS
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (SOSEQS-S, DSOSEQ-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     DSOSEQ solves a system of N simultaneous nonlinear equations.
C     See the comments in the interfacing routine DSOS for a more
C     detailed description of some of the items in the calling list.
C
C **********************************************************************
C   -Input-
C
C     FNC- Function subprogram which evaluates the equations
C     N  -number of equations
C     S  -Solution vector of initial guesses
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
C     IPRINT-Internal printing parameter. You must set IPRINT=-1 if you
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
C   -Output-
C     S    -Solution vector
C     IFLAG-Status indicator flag
C     MXIT-The actual number of iterations performed
C     FMAX-Residual norm
C     C   -Upper unit triangular matrix which approximates the
C          forward triangularization of the full Jacobian matrix.
C          Stored in a vector with dimension at least N*(N+1)/2.
C     B   -Contains the residuals (function values) divided
C          by the corresponding components of the P vector
C     P   -Array used to store the partial derivatives. After
C          each iteration P(K) contains the maximal derivative
C          occurring in the K-th reduced equation.
C     TEMP-Array used to store the previous solution iterate.
C     X   -Solution vector. Contains the values achieved on the
C          last iteration loop upon exit from DSOS.
C     Y   -Array containing the solution increments.
C     FAC -Array containing factors used in computing numerical
C          derivatives.
C     IS  -Records the pivotal information (column interchanges)
C
C **********************************************************************
C *** Three machine dependent parameters appear in this subroutine.
C
C *** The smallest positive magnitude, zero, is defined by the function
C *** routine D1MACH(1).
C
C *** URO, the computer unit roundoff value, is defined by D1MACH(3) for
C *** machines that round or D1MACH(4) for machines that truncate.
C *** URO is the smallest positive number such that 1.+URO  .GT.  1.
C
C *** The output tape unit number, LOUN, is defined by the function
C *** I1MACH(2).
C **********************************************************************
C
C***SEE ALSO  DSOS
C***ROUTINES CALLED  D1MACH, DSOSSL, I1MACH
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DSOSEQ
C
C
      INTEGER I1MACH
      DOUBLE PRECISION D1MACH
      INTEGER IC, ICR, IFLAG, IPRINT, IS(*), ISJ, ISV, IT, ITEM, ITRY,
     1     J, JK, JS, K, KD, KJ, KK, KM1, KN, KSV, L, LOUN, LS, M, MIT,
     2     MM, MXIT, N, NC, NCJS, NP1, NSRI, NSRRC
      DOUBLE PRECISION ATOLX, B(*), C(*), CSV, F, FAC(*), FACT, FDIF,
     1     FMAX, FMIN, FMXS, FN1, FN2, FNC, FP, H, HX, P(*), PMAX, RE,
     2     RTOLX, S(*), SRURO, TEMP(*), TEST, TOLF, URO, X(*), XNORM,
     3     Y(*), YJ, YN1, YN2, YN3, YNORM, YNS, ZERO
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 430
C        BEGIN BLOCK PERMITTING ...EXITS TO 410
C           BEGIN BLOCK PERMITTING ...EXITS TO 390
C***FIRST EXECUTABLE STATEMENT  DSOSEQ
               URO = D1MACH(4)
               LOUN = I1MACH(2)
               ZERO = D1MACH(1)
               RE = MAX(RTOLX,URO)
               SRURO = SQRT(URO)
C
               IFLAG = 0
               NP1 = N + 1
               ICR = 0
               IC = 0
               ITRY = NCJS
               YN1 = 0.0D0
               YN2 = 0.0D0
               YN3 = 0.0D0
               YNS = 0.0D0
               MIT = 0
               FN1 = 0.0D0
               FN2 = 0.0D0
               FMXS = 0.0D0
C
C              INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND
C              SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE.
C
               DO 10 K = 1, N
                  IS(K) = K
                  X(K) = S(K)
                  TEMP(K) = X(K)
   10          CONTINUE
C
C
C              *********************************************************
C              **** BEGIN PRINCIPAL ITERATION LOOP  ****
C              *********************************************************
C
               DO 380 M = 1, MXIT
C                 BEGIN BLOCK PERMITTING ...EXITS TO 350
C                    BEGIN BLOCK PERMITTING ...EXITS TO 240
C
                        DO 20 K = 1, N
                           FAC(K) = SRURO
   20                   CONTINUE
C
   30                   CONTINUE
C                          BEGIN BLOCK PERMITTING ...EXITS TO 180
                              KN = 1
                              FMAX = 0.0D0
C
C
C                             ******** BEGIN SUBITERATION LOOP DEFINING
C                             THE LINEARIZATION OF EACH ********
C                             EQUATION WHICH RESULTS IN THE CONSTRUCTION
C                             OF AN UPPER ******** TRIANGULAR MATRIX
C                             APPROXIMATING THE FORWARD ********
C                             TRIANGULARIZATION OF THE FULL JACOBIAN
C                             MATRIX
C
                              DO 170 K = 1, N
C                                BEGIN BLOCK PERMITTING ...EXITS TO 160
                                    KM1 = K - 1
C
C                                   BACK-SOLVE A TRIANGULAR LINEAR
C                                   SYSTEM OBTAINING IMPROVED SOLUTION
C                                   VALUES FOR K-1 OF THE VARIABLES FROM
C                                   THE FIRST K-1 EQUATIONS. THESE
C                                   VARIABLES ARE THEN ELIMINATED FROM
C                                   THE K-TH EQUATION.
C
                                    IF (KM1 .EQ. 0) GO TO 50
                                       CALL DSOSSL(K,N,KM1,Y,C,B,KN)
                                       DO 40 J = 1, KM1
                                          JS = IS(J)
                                          X(JS) = TEMP(JS) + Y(J)
   40                                  CONTINUE
   50                               CONTINUE
C
C
C                                   EVALUATE THE K-TH EQUATION AND THE
C                                   INTERMEDIATE COMPUTATION FOR THE MAX
C                                   NORM OF THE RESIDUAL VECTOR.
C
                                    F = FNC(X,K)
                                    FMAX = MAX(FMAX,ABS(F))
C
C                                   IF WE WISH TO PERFORM SEVERAL
C                                   ITERATIONS USING A FIXED
C                                   FACTORIZATION OF AN APPROXIMATE
C                                   JACOBIAN,WE NEED ONLY UPDATE THE
C                                   CONSTANT VECTOR.
C
C                                ...EXIT
                                    IF (ITRY .LT. NCJS) GO TO 160
C
C
                                    IT = 0
C
C                                   COMPUTE PARTIAL DERIVATIVES THAT ARE
C                                   REQUIRED IN THE LINEARIZATION OF THE
C                                   K-TH REDUCED EQUATION
C
                                    DO 90 J = K, N
                                       ITEM = IS(J)
                                       HX = X(ITEM)
                                       H = FAC(ITEM)*HX
                                       IF (ABS(H) .LE. ZERO)
     1                                    H = FAC(ITEM)
                                       X(ITEM) = HX + H
                                       IF (KM1 .EQ. 0) GO TO 70
                                          Y(J) = H
                                          CALL DSOSSL(K,N,J,Y,C,B,KN)
                                          DO 60 L = 1, KM1
                                             LS = IS(L)
                                             X(LS) = TEMP(LS) + Y(L)
   60                                     CONTINUE
   70                                  CONTINUE
                                       FP = FNC(X,K)
                                       X(ITEM) = HX
                                       FDIF = FP - F
                                       IF (ABS(FDIF) .GT. URO*ABS(F))
     1                                    GO TO 80
                                          FDIF = 0.0D0
                                          IT = IT + 1
   80                                  CONTINUE
                                       P(J) = FDIF/H
   90                               CONTINUE
C
                                    IF (IT .LE. (N - K)) GO TO 110
C
C                                      ALL COMPUTED PARTIAL DERIVATIVES
C                                      OF THE K-TH EQUATION ARE
C                                      EFFECTIVELY ZERO.TRY LARGER
C                                      PERTURBATIONS OF THE INDEPENDENT
C                                      VARIABLES.
C
                                       DO 100 J = K, N
                                          ISJ = IS(J)
                                          FACT = 100.0D0*FAC(ISJ)
C           ..............................EXIT
                                          IF (FACT .GT. 1.0D10)
     1                                       GO TO 390
                                          FAC(ISJ) = FACT
  100                                  CONTINUE
C                          ............EXIT
                                       GO TO 180
  110                               CONTINUE
C
C                                ...EXIT
                                    IF (K .EQ. N) GO TO 160
C
C                                   ACHIEVE A PIVOTING EFFECT BY
C                                   CHOOSING THE MAXIMAL DERIVATIVE
C                                   ELEMENT
C
                                    PMAX = 0.0D0
                                    DO 130 J = K, N
                                       TEST = ABS(P(J))
                                       IF (TEST .LE. PMAX) GO TO 120
                                          PMAX = TEST
                                          ISV = J
  120                                  CONTINUE
  130                               CONTINUE
C           ........................EXIT
                                    IF (PMAX .EQ. 0.0D0) GO TO 390
C
C                                   SET UP THE COEFFICIENTS FOR THE K-TH
C                                   ROW OF THE TRIANGULAR LINEAR SYSTEM
C                                   AND SAVE THE PARTIAL DERIVATIVE OF
C                                   LARGEST MAGNITUDE
C
                                    PMAX = P(ISV)
                                    KK = KN
                                    DO 140 J = K, N
                                       IF (J .NE. ISV)
     1                                    C(KK) = -P(J)/PMAX
                                       KK = KK + 1
  140                               CONTINUE
                                    P(K) = PMAX
C
C
C                                ...EXIT
                                    IF (ISV .EQ. K) GO TO 160
C
C                                   INTERCHANGE THE TWO COLUMNS OF C
C                                   DETERMINED BY THE PIVOTAL STRATEGY
C
                                    KSV = IS(K)
                                    IS(K) = IS(ISV)
                                    IS(ISV) = KSV
C
                                    KD = ISV - K
                                    KJ = K
                                    DO 150 J = 1, K
                                       CSV = C(KJ)
                                       JK = KJ + KD
                                       C(KJ) = C(JK)
                                       C(JK) = CSV
                                       KJ = KJ + N - J
  150                               CONTINUE
  160                            CONTINUE
C
                                 KN = KN + NP1 - K
C
C                                STORE THE COMPONENTS FOR THE CONSTANT
C                                VECTOR
C
                                 B(K) = -F/P(K)
C
  170                         CONTINUE
C                       ......EXIT
                              GO TO 190
  180                      CONTINUE
                        GO TO 30
  190                   CONTINUE
C
C                       ********
C                       ******** END OF LOOP CREATING THE TRIANGULAR
C                       LINEARIZATION MATRIX
C                       ********
C
C
C                        SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW
C                        SOLUTION APPROXIMATION AND OBTAIN THE SOLUTION
C                        INCREMENT NORM.
C
                        KN = KN - 1
                        Y(N) = B(N)
                        IF (N .GT. 1) CALL DSOSSL(N,N,N,Y,C,B,KN)
                        XNORM = 0.0D0
                        YNORM = 0.0D0
                        DO 200 J = 1, N
                           YJ = Y(J)
                           YNORM = MAX(YNORM,ABS(YJ))
                           JS = IS(J)
                           X(JS) = TEMP(JS) + YJ
                           XNORM = MAX(XNORM,ABS(X(JS)))
  200                   CONTINUE
C
C
C                       PRINT INTERMEDIATE SOLUTION ITERATES AND
C                       RESIDUAL NORM IF DESIRED
C
                        IF (IPRINT .NE. (-1)) GO TO 220
                           MM = M - 1
                           WRITE (LOUN,210) FMAX,MM,(X(J), J = 1, N)
  210                      FORMAT ('0RESIDUAL NORM =', D9.2, / 1X,
     1                             'SOLUTION ITERATE (', I3, ')', /
     2                             (1X, 5D26.14))
  220                   CONTINUE
C
C                       TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE
C                       AND/OR ABSOLUTE ERROR COMPARISON ON SUCCESSIVE
C                       APPROXIMATIONS OF EACH SOLUTION VARIABLE)
C
                        DO 230 J = 1, N
                           JS = IS(J)
C                    ......EXIT
                           IF (ABS(Y(J)) .GT. RE*ABS(X(JS)) + ATOLX)
     1                        GO TO 240
  230                   CONTINUE
                        IF (FMAX .LE. FMXS) IFLAG = 1
  240                CONTINUE
C
C                    TEST FOR CONVERGENCE TO A SOLUTION BASED ON
C                    RESIDUALS
C
                     IF (FMAX .LE. TOLF) IFLAG = IFLAG + 2
C        ............EXIT
                     IF (IFLAG .GT. 0) GO TO 410
C
C
                     IF (M .GT. 1) GO TO 250
                        FMIN = FMAX
                     GO TO 330
  250                CONTINUE
C                       BEGIN BLOCK PERMITTING ...EXITS TO 320
C
C                          SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM.
C
                           IF (FMAX .GE. FMIN) GO TO 270
                              MIT = M + 1
                              YN1 = YNORM
                              YN2 = YNS
                              FN1 = FMXS
                              FMIN = FMAX
                              DO 260 J = 1, N
                                 S(J) = X(J)
  260                         CONTINUE
                              IC = 0
  270                      CONTINUE
C
C                          TEST FOR LIMITING PRECISION CONVERGENCE. VERY
C                          SLOWLY CONVERGENT PROBLEMS MAY ALSO BE
C                          DETECTED.
C
                           IF (YNORM .GT. SRURO*XNORM) GO TO 290
                           IF (FMAX .LT. 0.2D0*FMXS
     1                         .OR. FMAX .GT. 5.0D0*FMXS) GO TO 290
                           IF (YNORM .LT. 0.2D0*YNS
     1                         .OR. YNORM .GT. 5.0D0*YNS) GO TO 290
                              ICR = ICR + 1
                              IF (ICR .GE. NSRRC) GO TO 280
                                 IC = 0
C                       .........EXIT
                                 GO TO 320
  280                         CONTINUE
                              IFLAG = 4
                              FMAX = FMIN
C     ........................EXIT
                              GO TO 430
  290                      CONTINUE
                           ICR = 0
C
C                          TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME.
C
                           IF (YNORM .GT. 2.0D0*YNS
     1                         .OR. FMAX .GT. 2.0D0*FMXS) GO TO 300
                              IC = 0
                           GO TO 310
  300                      CONTINUE
                              IC = IC + 1
C                       ......EXIT
                              IF (IC .LT. NSRI) GO TO 320
                              IFLAG = 7
C        .....................EXIT
                              GO TO 410
  310                      CONTINUE
  320                   CONTINUE
  330                CONTINUE
C
C                    CHECK TO SEE IF NEXT ITERATION CAN USE THE OLD
C                    JACOBIAN FACTORIZATION
C
                     ITRY = ITRY - 1
                     IF (ITRY .EQ. 0) GO TO 340
                     IF (20.0D0*YNORM .GT. XNORM) GO TO 340
                     IF (YNORM .GT. 2.0D0*YNS) GO TO 340
C                 ......EXIT
                        IF (FMAX .LT. 2.0D0*FMXS) GO TO 350
  340                CONTINUE
                     ITRY = NCJS
  350             CONTINUE
C
C                 SAVE THE CURRENT SOLUTION APPROXIMATION AND THE
C                 RESIDUAL AND SOLUTION INCREMENT NORMS FOR USE IN THE
C                 NEXT ITERATION.
C
                  DO 360 J = 1, N
                     TEMP(J) = X(J)
  360             CONTINUE
                  IF (M .NE. MIT) GO TO 370
                     FN2 = FMAX
                     YN3 = YNORM
  370             CONTINUE
                  FMXS = FMAX
                  YNS = YNORM
C
C
  380          CONTINUE
C
C              *********************************************************
C              **** END OF PRINCIPAL ITERATION LOOP ****
C              *********************************************************
C
C
C               TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED.
               M = MXIT
               IFLAG = 5
               IF (YN1 .GT. 10.0D0*YN2 .OR. YN3 .GT. 10.0D0*YN1)
     1            IFLAG = 6
               IF (FN1 .GT. 5.0D0*FMIN .OR. FN2 .GT. 5.0D0*FMIN)
     1            IFLAG = 6
               IF (FMAX .GT. 5.0D0*FMIN) IFLAG = 6
C        ......EXIT
               GO TO 410
  390       CONTINUE
C
C
C           A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR.
            IFLAG = 8
            DO 400 J = 1, N
               S(J) = TEMP(J)
  400       CONTINUE
C     ......EXIT
            GO TO 430
  410    CONTINUE
C
C
         DO 420 J = 1, N
            S(J) = X(J)
  420    CONTINUE
  430 CONTINUE
C
C
      MXIT = M
      RETURN
      END
