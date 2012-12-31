*DECK DNSQ
      SUBROUTINE DNSQ (FCN, JAC, IOPT, N, X, FVEC, FJAC, LDFJAC, XTOL,
     +   MAXFEV, ML, MU, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, NFEV,
     +   NJEV, R, LR, QTF, WA1, WA2, WA3, WA4)
C***BEGIN PROLOGUE  DNSQ
C***PURPOSE  Find a zero of a system of a N nonlinear functions in N
C            variables by a modification of the Powell hybrid method.
C***LIBRARY   SLATEC
C***CATEGORY  F2A
C***TYPE      DOUBLE PRECISION (SNSQ-S, DNSQ-D)
C***KEYWORDS  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS
C***AUTHOR  Hiebert, K. L. (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of DNSQ is to find a zero of a system of N nonlinear
C       functions in N variables by a modification of the Powell
C       hybrid method.  The user must provide a subroutine which
C       calculates the functions.  The user has the option of either to
C       provide a subroutine which calculates the Jacobian or to let the
C       code calculate it by a forward-difference approximation.
C       This code is the combination of the MINPACK codes (Argonne)
C       HYBRD and HYBRDJ.
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
C      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
C      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
C       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
C       DOUBLE PRECISION XTOL,EPSFCN,FACTOR
C       DOUBLE PRECISION
C       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
C      *     WA1(N),WA2(N),WA3(N),WA4(N)
C       EXTERNAL FCN,JAC
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to DNSQ and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from DNSQ.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  FCN must be declared in an EXTERNAL statement
C         in the user calling program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of DNSQ.  In this case set
C         IFLAG to a negative integer.
C
C       JAC is the name of the user-supplied subroutine which calculates
C         the Jacobian.  If IOPT=1, then JAC must be declared in an
C         EXTERNAL statement in the user calling program, and should be
C         written as follows.
C
C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C         INTEGER N,LDFJAC,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
C         ----------
C         Calculate the Jacobian at X and return this
C         matrix in FJAC.  FVEC contains the function
C         values at X and should not be altered.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by JAC unless the
C         user wants to terminate execution of DNSQ.  In this case set
C         IFLAG to a negative integer.
C
C         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=1, then the user must supply the
C         Jacobian through the subroutine JAC.  If IOPT=2, then the
C         code will approximate the Jacobian by forward-differencing.
C
C       N is a positive integer input variable set to the number of
C         functions and variables.
C
C       X is an array of length N.  On input X must contain an initial
C         estimate of the solution vector.  On output X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length N which contains the functions
C         evaluated at the output X.
C
C       FJAC is an output N by N array which contains the orthogonal
C         matrix Q produced by the QR factorization of the final
C         approximate Jacobian.
C
C       LDFJAC is a positive integer input variable not less than N
C         which specifies the leading dimension of the array FJAC.
C
C       XTOL is a nonnegative input variable.  Termination occurs when
C         the relative error between two consecutive iterates is at most
C         XTOL.  Therefore, XTOL measures the relative error desired in
C         the approximate solution.  Section 4 contains more details
C         about XTOL.
C
C       MAXFEV is a positive integer input variable.  Termination occurs
C         when the number of calls to FCN is at least MAXFEV by the end
C         of an iteration.
C
C       ML is a nonnegative integer input variable which specifies the
C         number of subdiagonals within the band of the Jacobian matrix.
C         If the Jacobian is not banded or IOPT=1, set ML to at
C         least N - 1.
C
C       MU is a nonnegative integer input variable which specifies the
C         number of superdiagonals within the band of the Jacobian
C         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
C         least N - 1.
C
C       EPSFCN is an input variable used in determining a suitable step
C         for the forward-difference approximation.  This approximation
C         assumes that the relative errors in the functions are of the
C         order of EPSFCN.  If EPSFCN is less than the machine
C         precision, it is assumed that the relative errors in the
C         functions are of the order of the machine precision.  If
C         IOPT=1, then EPSFCN can be ignored (treat it as a dummy
C         argument).
C
C       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
C         internally set.  If MODE = 2, DIAG must contain positive
C         entries that serve as implicit (multiplicative) scale factors
C         for the variables.
C
C       MODE is an integer input variable.  If MODE = 1, the variables
C         will be scaled internally.  If MODE = 2, the scaling is
C         specified by the input DIAG.  Other values of MODE are
C         equivalent to MODE = 1.
C
C       FACTOR is a positive input variable used in determining the
C         initial step bound.  This bound is set to the product of
C         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to
C         FACTOR itself.  In most cases FACTOR should lie in the
C         interval (.1,100.).  100. is a generally recommended value.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. appropriate
C         print statements must be added to FCN(see example).  If NPRINT
C         is not positive, no special calls of FCN with IFLAG = 0 are
C         made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  Improper input parameters.
C
C         INFO = 1  Relative error between two consecutive iterates is
C                   at most XTOL.
C
C         INFO = 2  Number of calls to FCN has reached or exceeded
C                   MAXFEV.
C
C         INFO = 3  XTOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 4  Iteration is not making good progress, as measured
C                   by the improvement from the last five Jacobian
C                   evaluations.
C
C         INFO = 5  Iteration is not making good progress, as measured
C                   by the improvement from the last ten iterations.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       NFEV is an integer output variable set to the number of calls to
C         FCN.
C
C       NJEV is an integer output variable set to the number of calls to
C         JAC. (If IOPT=2, then NJEV is set to zero.)
C
C       R is an output array of length LR which contains the upper
C         triangular matrix produced by the QR factorization of the
C         final approximate Jacobian, stored rowwise.
C
C       LR is a positive integer input variable not less than
C         (N*(N+1))/2.
C
C       QTF is an output array of length N which contains the vector
C         (Q transpose)*FVEC.
C
C       WA1, WA2, WA3, and WA4 are work arrays of length N.
C
C
C 4. Successful completion.
C
C       The accuracy of DNSQ is controlled by the convergence parameter
C       XTOL.  This parameter is used in a test which makes a comparison
C       between the approximation X and a solution XSOL.  DNSQ
C       terminates when the test is satisfied.  If the convergence
C       parameter is less than the machine precision (as defined by the
C       function D1MACH(4)), then DNSQ only attempts to satisfy the test
C       defined by the machine precision.  Further progress is not
C       usually possible.
C
C       The test assumes that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently.  If these conditions
C       are not satisfied, then DNSQ may incorrectly indicate
C       convergence.  The coding of the Jacobian can be checked by the
C       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2,
C       then the validity of the answer can be checked, for example, by
C       rerunning DNSQ with a tighter tolerance.
C
C       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
C         vector Z and D is the diagonal matrix whose entries are
C         defined by the array DIAG, then this test attempts to
C         guarantee that
C
C               DENORM(D*(X-XSOL)) .LE. XTOL*DENORM(D*XSOL).
C
C         If this condition is satisfied with XTOL = 10**(-K), then the
C         larger components of D*X have K significant decimal digits and
C         INFO is set to 1.  There is a danger that the smaller
C         components of D*X may have large relative errors, but the fast
C         rate of convergence of DNSQ usually avoids this possibility.
C         Unless high precision solutions are required, the recommended
C         value for XTOL is the square root of the machine precision.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of DNSQ can be due to improper input
C       parameters, arithmetic interrupts, an excessive number of
C       function evaluations, or lack of good progress.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1,
C         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or
C         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0,
C         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by DNSQ.  In this
C         case, it may be possible to remedy the situation by rerunning
C         DNSQ with a smaller value of FACTOR.
C
C       Excessive Number of Function Evaluations.  A reasonable value
C         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
C         If the number of calls to FCN reaches MAXFEV, then this
C         indicates that the routine is converging very slowly as
C         measured by the progress of FVEC, and INFO is set to 2. This
C         situation should be unusual because, as indicated below, lack
C         of good progress is usually diagnosed earlier by DNSQ,
C         causing termination with info = 4 or INFO = 5.
C
C       Lack of Good Progress.  DNSQ searches for a zero of the system
C         by minimizing the sum of the squares of the functions.  In so
C         doing, it can become trapped in a region where the minimum
C         does not correspond to a zero of the system and, in this
C         situation, the iteration eventually fails to make good
C         progress.  In particular, this will happen if the system does
C         not have a zero.  If the system has a zero, rerunning DNSQ
C         from a different starting point may be helpful.
C
C
C 6. Characteristics of The Algorithm.
C
C       DNSQ is a modification of the Powell Hybrid method.  Two of its
C       main characteristics involve the choice of the correction as a
C       convex combination of the Newton and scaled gradient directions,
C       and the updating of the Jacobian by the rank-1 method of
C       Broyden.  The choice of the correction guarantees (under
C       reasonable conditions) global convergence for starting points
C       far from the solution and a fast rate of convergence.  The
C       Jacobian is calculated at the starting point by either the
C       user-supplied subroutine or a forward-difference approximation,
C       but it is not recalculated until the rank-1 method fails to
C       produce satisfactory progress.
C
C       Timing.  The time required by DNSQ to solve a given problem
C         depends on N, the behavior of the functions, the accuracy
C         requested, and the starting point.  The number of arithmetic
C         operations needed by DNSQ is about 11.5*(N**2) to process
C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
C         to process each evaluation of the Jacobian (call to JAC,
C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
C         the timing of DNSQ will be strongly influenced by the time
C         spent in FCN and JAC.
C
C       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision
C         storage locations, in addition to the storage required by the
C         program.  There are no internally declared storage arrays.
C
C *Long Description:
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), ..., X(9),
C       which solve the system of tridiagonal equations
C
C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
C                                   -X(8) + (3-2*X(9))*X(9) = -1
C C     **********
C
C       PROGRAM TEST
C C
C C     Driver for DNSQ example.
C C
C       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
C      *        NWRITE
C       DOUBLE PRECISION XTOL,EPSFCN,FACTOR,FNORM
C       DOUBLE PRECISION X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
C      *     WA1(9),WA2(9),WA3(9),WA4(9)
C       DOUBLE PRECISION DENORM,D1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 2
C       N = 9
C C
C C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
C C
C       DO 10 J = 1, 9
C          X(J) = -1.E0
C    10    CONTINUE
C C
C       LDFJAC = 9
C       LR = 45
C C
C C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C C     THIS IS THE RECOMMENDED SETTING.
C C
C       XTOL = SQRT(D1MACH(4))
C C
C       MAXFEV = 2000
C       ML = 1
C       MU = 1
C       EPSFCN = 0.E0
C       MODE = 2
C       DO 20 J = 1, 9
C          DIAG(J) = 1.E0
C    20    CONTINUE
C       FACTOR = 1.E2
C       NPRINT = 0
C C
C       CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
C      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
C      *           R,LR,QTF,WA1,WA2,WA3,WA4)
C       FNORM = DENORM(N,FVEC)
C       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
C       END
C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
C       INTEGER N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(N)
C       INTEGER K
C       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.
C C
C       RETURN
C     5 CONTINUE
C       DO 10 K = 1, N
C          TEMP = (THREE - TWO*X(K))*X(K)
C          TEMP1 = ZERO
C          IF (K .NE. 1) TEMP1 = X(K-1)
C          TEMP2 = ZERO
C          IF (K .NE. N) TEMP2 = X(K+1)
C          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
C    10    CONTINUE
C       RETURN
C       END
C
C       Results obtained with different compilers or machines
C       may be slightly different.
C
C       Final L2 norm of the residuals  0.1192636E-07
C
C       Number of function evaluations        14
C
C       Exit parameter                         1
C
C       Final approximate solution
C
C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
C
C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
C                 tions. In Numerical Methods for Nonlinear Algebraic
C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
C                 1988.
C***ROUTINES CALLED  D1MACH, D1MPYQ, D1UPDT, DDOGLG, DENORM, DFDJC1,
C                    DQFORM, DQRFAC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNSQ
      DOUBLE PRECISION D1MACH,DENORM
      INTEGER I, IFLAG, INFO, IOPT, ITER, IWA(1), J, JM1, L, LDFJAC,
     1     LR, MAXFEV, ML, MODE, MU, N, NCFAIL, NCSUC, NFEV, NJEV,
     2     NPRINT, NSLOW1, NSLOW2
      DOUBLE PRECISION ACTRED, DELTA, DIAG(*), EPSFCN, EPSMCH, FACTOR,
     1     FJAC(LDFJAC,*), FNORM, FNORM1, FVEC(*), ONE, P0001, P001,
     2     P1, P5, PNORM, PRERED, QTF(*), R(*), RATIO, SUM, TEMP,
     3     WA1(*), WA2(*), WA3(*), WA4(*), X(*), XNORM, XTOL, ZERO
      EXTERNAL FCN
      LOGICAL JEVAL,SING
      SAVE ONE, P1, P5, P001, P0001, ZERO
      DATA ONE,P1,P5,P001,P0001,ZERO
     1     /1.0D0,1.0D-1,5.0D-1,1.0D-3,1.0D-4,0.0D0/
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 320
C***FIRST EXECUTABLE STATEMENT  DNSQ
         EPSMCH = D1MACH(4)
C
         INFO = 0
         IFLAG = 0
         NFEV = 0
         NJEV = 0
C
C        CHECK THE INPUT PARAMETERS FOR ERRORS.
C
C     ...EXIT
         IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0
     1       .OR. XTOL .LT. ZERO .OR. MAXFEV .LE. 0 .OR. ML .LT. 0
     2       .OR. MU .LT. 0 .OR. FACTOR .LE. ZERO .OR. LDFJAC .LT. N
     3       .OR. LR .LT. (N*(N + 1))/2) GO TO 320
         IF (MODE .NE. 2) GO TO 20
            DO 10 J = 1, N
C     .........EXIT
               IF (DIAG(J) .LE. ZERO) GO TO 320
   10       CONTINUE
   20    CONTINUE
C
C        EVALUATE THE FUNCTION AT THE STARTING POINT
C        AND CALCULATE ITS NORM.
C
         IFLAG = 1
         CALL FCN(N,X,FVEC,IFLAG)
         NFEV = 1
C     ...EXIT
         IF (IFLAG .LT. 0) GO TO 320
         FNORM = DENORM(N,FVEC)
C
C        INITIALIZE ITERATION COUNTER AND MONITORS.
C
         ITER = 1
         NCSUC = 0
         NCFAIL = 0
         NSLOW1 = 0
         NSLOW2 = 0
C
C        BEGINNING OF THE OUTER LOOP.
C
   30    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 90
               JEVAL = .TRUE.
C
C              CALCULATE THE JACOBIAN MATRIX.
C
               IF (IOPT .EQ. 2) GO TO 40
C
C                 USER SUPPLIES JACOBIAN
C
                  CALL JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
                  NJEV = NJEV + 1
               GO TO 50
   40          CONTINUE
C
C                 CODE APPROXIMATES THE JACOBIAN
C
                  IFLAG = 2
                  CALL DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,
     1                        EPSFCN,WA1,WA2)
                  NFEV = NFEV + MIN(ML+MU+1,N)
   50          CONTINUE
C
C     .........EXIT
               IF (IFLAG .LT. 0) GO TO 320
C
C              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
               CALL DQRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
C
C              ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
C           ...EXIT
               IF (ITER .NE. 1) GO TO 90
               IF (MODE .EQ. 2) GO TO 70
                  DO 60 J = 1, N
                     DIAG(J) = WA2(J)
                     IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   60             CONTINUE
   70          CONTINUE
C
C              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED
C              X AND INITIALIZE THE STEP BOUND DELTA.
C
               DO 80 J = 1, N
                  WA3(J) = DIAG(J)*X(J)
   80          CONTINUE
               XNORM = DENORM(N,WA3)
               DELTA = FACTOR*XNORM
               IF (DELTA .EQ. ZERO) DELTA = FACTOR
   90       CONTINUE
C
C           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
C
            DO 100 I = 1, N
               QTF(I) = FVEC(I)
  100       CONTINUE
            DO 140 J = 1, N
               IF (FJAC(J,J) .EQ. ZERO) GO TO 130
                  SUM = ZERO
                  DO 110 I = J, N
                     SUM = SUM + FJAC(I,J)*QTF(I)
  110             CONTINUE
                  TEMP = -SUM/FJAC(J,J)
                  DO 120 I = J, N
                     QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
C
C           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
C
            SING = .FALSE.
            DO 170 J = 1, N
               L = J
               JM1 = J - 1
               IF (JM1 .LT. 1) GO TO 160
               DO 150 I = 1, JM1
                  R(L) = FJAC(I,J)
                  L = L + N - I
  150          CONTINUE
  160          CONTINUE
               R(L) = WA1(J)
               IF (WA1(J) .EQ. ZERO) SING = .TRUE.
  170       CONTINUE
C
C           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
C
            CALL DQFORM(N,N,FJAC,LDFJAC,WA1)
C
C           RESCALE IF NECESSARY.
C
            IF (MODE .EQ. 2) GO TO 190
               DO 180 J = 1, N
                  DIAG(J) = MAX(DIAG(J),WA2(J))
  180          CONTINUE
  190       CONTINUE
C
C           BEGINNING OF THE INNER LOOP.
C
  200       CONTINUE
C
C              IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
               IF (NPRINT .LE. 0) GO TO 210
                  IFLAG = 0
                  IF (MOD(ITER-1,NPRINT) .EQ. 0)
     1               CALL FCN(N,X,FVEC,IFLAG)
C     ............EXIT
                  IF (IFLAG .LT. 0) GO TO 320
  210          CONTINUE
C
C              DETERMINE THE DIRECTION P.
C
               CALL DDOGLG(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C
C              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
               DO 220 J = 1, N
                  WA1(J) = -WA1(J)
                  WA2(J) = X(J) + WA1(J)
                  WA3(J) = DIAG(J)*WA1(J)
  220          CONTINUE
               PNORM = DENORM(N,WA3)
C
C              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
               IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)
C
C              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
               IFLAG = 1
               CALL FCN(N,WA2,WA4,IFLAG)
               NFEV = NFEV + 1
C     .........EXIT
               IF (IFLAG .LT. 0) GO TO 320
               FNORM1 = DENORM(N,WA4)
C
C              COMPUTE THE SCALED ACTUAL REDUCTION.
C
               ACTRED = -ONE
               IF (FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C              COMPUTE THE SCALED PREDICTED REDUCTION.
C
               L = 1
               DO 240 I = 1, N
                  SUM = ZERO
                  DO 230 J = I, N
                     SUM = SUM + R(L)*WA1(J)
                     L = L + 1
  230             CONTINUE
                  WA3(I) = QTF(I) + SUM
  240          CONTINUE
               TEMP = DENORM(N,WA3)
               PRERED = ZERO
               IF (TEMP .LT. FNORM) PRERED = ONE - (TEMP/FNORM)**2
C
C              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C              REDUCTION.
C
               RATIO = ZERO
               IF (PRERED .GT. ZERO) RATIO = ACTRED/PRERED
C
C              UPDATE THE STEP BOUND.
C
               IF (RATIO .GE. P1) GO TO 250
                  NCSUC = 0
                  NCFAIL = NCFAIL + 1
                  DELTA = P5*DELTA
               GO TO 260
  250          CONTINUE
                  NCFAIL = 0
                  NCSUC = NCSUC + 1
                  IF (RATIO .GE. P5 .OR. NCSUC .GT. 1)
     1               DELTA = MAX(DELTA,PNORM/P5)
                  IF (ABS(RATIO-ONE) .LE. P1) DELTA = PNORM/P5
  260          CONTINUE
C
C              TEST FOR SUCCESSFUL ITERATION.
C
               IF (RATIO .LT. P0001) GO TO 280
C
C                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
                  DO 270 J = 1, N
                     X(J) = WA2(J)
                     WA2(J) = DIAG(J)*X(J)
                     FVEC(J) = WA4(J)
  270             CONTINUE
                  XNORM = DENORM(N,WA2)
                  FNORM = FNORM1
                  ITER = ITER + 1
  280          CONTINUE
C
C              DETERMINE THE PROGRESS OF THE ITERATION.
C
               NSLOW1 = NSLOW1 + 1
               IF (ACTRED .GE. P001) NSLOW1 = 0
               IF (JEVAL) NSLOW2 = NSLOW2 + 1
               IF (ACTRED .GE. P1) NSLOW2 = 0
C
C              TEST FOR CONVERGENCE.
C
               IF (DELTA .LE. XTOL*XNORM .OR. FNORM .EQ. ZERO) INFO = 1
C     .........EXIT
               IF (INFO .NE. 0) GO TO 320
C
C              TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
               IF (NFEV .GE. MAXFEV) INFO = 2
               IF (P1*MAX(P1*DELTA,PNORM) .LE. EPSMCH*XNORM) INFO = 3
               IF (NSLOW2 .EQ. 5) INFO = 4
               IF (NSLOW1 .EQ. 10) INFO = 5
C     .........EXIT
               IF (INFO .NE. 0) GO TO 320
C
C              CRITERION FOR RECALCULATING JACOBIAN
C
C           ...EXIT
               IF (NCFAIL .EQ. 2) GO TO 310
C
C              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
C              AND UPDATE QTF IF NECESSARY.
C
               DO 300 J = 1, N
                  SUM = ZERO
                  DO 290 I = 1, N
                     SUM = SUM + FJAC(I,J)*WA4(I)
  290             CONTINUE
                  WA2(J) = (SUM - WA3(J))/PNORM
                  WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
                  IF (RATIO .GE. P0001) QTF(J) = SUM
  300          CONTINUE
C
C              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
C
               CALL D1UPDT(N,N,R,LR,WA1,WA2,WA3,SING)
               CALL D1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)
               CALL D1MPYQ(1,N,QTF,1,WA2,WA3)
C
C              END OF THE INNER LOOP.
C
               JEVAL = .FALSE.
            GO TO 200
  310       CONTINUE
C
C           END OF THE OUTER LOOP.
C
         GO TO 30
  320 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(N,X,FVEC,IFLAG)
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 2) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
      IF (INFO .EQ. 3) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'XTOL TOO SMALL. NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
      IF (INFO .GT. 4) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'ITERATION NOT MAKING GOOD PROGRESS.', 1, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE DNSQ.
C
      END
