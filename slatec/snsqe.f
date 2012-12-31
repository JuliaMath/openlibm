*DECK SNSQE
      SUBROUTINE SNSQE (FCN, JAC, IOPT, N, X, FVEC, TOL, NPRINT, INFO,
     +   WA, LWA)
C***BEGIN PROLOGUE  SNSQE
C***PURPOSE  An easy-to-use code to find a zero of a system of N
C            nonlinear functions in N variables by a modification of
C            the Powell hybrid method.
C***LIBRARY   SLATEC
C***CATEGORY  F2A
C***TYPE      SINGLE PRECISION (SNSQE-S, DNSQE-D)
C***KEYWORDS  EASY-TO-USE, NONLINEAR SQUARE SYSTEM,
C             POWELL HYBRID METHOD, ZEROS
C***AUTHOR  Hiebert, K. L., (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C
C       The purpose of SNSQE is to find a zero of a system of N non-
C       linear functions in N variables by a modification of the Powell
C       hybrid method.  This is done by using the more general nonlinear
C       equation solver SNSQ.  The user must provide a subroutine which
C       calculates the functions.  The user has the option of either to
C       provide a subroutine which calculates the Jacobian or to let the
C       code calculate it by a forward-difference approximation.  This
C       code is the combination of the MINPACK codes (Argonne) HYBRD1
C       and HYBRJ1.
C
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE SNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,
C      *                  WA,LWA)
C       INTEGER IOPT,N,NPRINT,INFO,LWA
C       REAL TOL
C       REAL X(N),FVEC(N),WA(LWA)
C       EXTERNAL FCN,JAC
C
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to SNSQE and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from SNSQE.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  FCN must be declared in an EXTERNAL statement
C         in the user calling program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         REAL X(N),FVEC(N)
C         ----------
C         Calculate the functions at X and
C         return this vector in FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of SNSQE.  In this case, set
C         IFLAG to a negative integer.
C
C       JAC is the name of the user-supplied subroutine which calculates
C         the Jacobian.  If IOPT=1, then JAC must be declared in an
C         EXTERNAL statement in the user calling program, and should be
C         written as follows.
C
C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C         INTEGER N,LDFJAC,IFLAG
C         REAL X(N),FVEC(N),FJAC(LDFJAC,N)
C         ----------
C         Calculate the Jacobian at X and return this
C         matrix in FJAC.  FVEC contains the function
C         values at X and should not be altered.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by JAC unless the
C         user wants to terminate execution of SNSQE.  In this case, set
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
C       X is an array of length N.  On input, X must contain an initial
C         estimate of the solution vector.  On output, X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length N which contains the functions
C         evaluated at the output X.
C
C       TOL is a non-negative input variable.  Termination occurs when
C         the algorithm estimates that the relative error between X and
C         the solution is at most TOL.  Section 4 contains more details
C         about TOL.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iteration thereafter and immediately prior
C         to return, with X and FVEC available for printing. Appropriate
C         print statements must be added to FCN (see example). If NPRINT
C         is not positive, no special calls of FCN with IFLAG = 0 are
C         made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  improper input parameters.
C
C         INFO = 1  algorithm estimates that the relative error between
C                   X and the solution is at most TOL.
C
C         INFO = 2  number of calls to FCN has reached or exceeded
C                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2.
C
C         INFO = 3  TOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 4  iteration is not making good progress.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       WA is a work array of length LWA.
C
C       LWA is a positive integer input variable not less than
C         (3*N**2+13*N))/2.
C
C
C 4. Successful Completion.
C
C       The accuracy of SNSQE is controlled by the convergence parame-
C       ter TOL.  This parameter is used in a test which makes a compar-
C       ison between the approximation X and a solution XSOL.  SNSQE
C       terminates when the test is satisfied.  If TOL is less than the
C       machine precision (as defined by the function R1MACH(4)), then
C       SNSQE attempts only to satisfy the test defined by the machine
C       precision.  Further progress is not usually possible.  Unless
C       high precision solutions are required, the recommended value
C       for TOL is the square root of the machine precision.
C
C       The test assumes that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian  coded consistently.  If these conditions
C       are not satisfied, SNSQE may incorrectly indicate convergence.
C       The coding of the Jacobian can be checked by the subroutine
C       CHKDER.  If the Jacobian is coded correctly or IOPT=2, then
C       the validity of the answer can be checked, for example, by
C       rerunning SNSQE with a tighter tolerance.
C
C       Convergence Test.  If ENORM(Z) denotes the Euclidean norm of a
C         vector Z, then this test attempts to guarantee that
C
C               ENORM(X-XSOL) .LE.  TOL*ENORM(XSOL).
C
C         If this condition is satisfied with TOL = 10**(-K), then the
C         larger components of X have K significant decimal digits and
C         INFO is set to 1.  There is a danger that the smaller compo-
C         nents of X may have large relative errors, but the fast rate
C         of convergence of SNSQE usually avoids this possibility.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of SNSQE can be due to improper input
C       parameters, arithmetic interrupts, an excessive number of func-
C       tion evaluations, errors in the functions, or lack of good prog-
C       ress.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, or
C         IOPT .GT. 2, or N .LE. 0, or TOL .LT. 0.E0, or
C         LWA .LT. (3*N**2+13*N)/2.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by SNSQE.  In this
C         case, it may be possible to remedy the situation by not evalu-
C         ating the functions here, but instead setting the components
C         of FVEC to numbers that exceed those in the initial FVEC.
C
C       Excessive Number of Function Evaluations.  If the number of
C         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for
C         IOPT=2, then this indicates that the routine is converging
C         very slowly as measured by the progress of FVEC, and INFO is
C         set to 2.  This situation should be unusual because, as
C         indicated below, lack of good progress is usually diagnosed
C         earlier by SNSQE, causing termination with INFO = 4.
C
C       Errors in the Functions.  When IOPT=2, the choice of step length
C         in the forward-difference approximation to the Jacobian
C         assumes that the relative errors in the functions are of the
C         order of the machine precision.  If this is not the case,
C         SNSQE may fail (usually with INFO = 4).  The user should
C         then either use SNSQ and set the step length or use IOPT=1
C         and supply the Jacobian.
C
C       Lack of Good Progress.  SNSQE searches for a zero of the system
C         by minimizing the sum of the squares of the functions.  In so
C         doing, it can become trapped in a region where the minimum
C         does not correspond to a zero of the system and, in this situ-
C         ation, the iteration eventually fails to make good progress.
C         In particular, this will happen if the system does not have a
C         zero.  If the system has a zero, rerunning SNSQE from a dif-
C         ferent starting point may be helpful.
C
C
C 6. Characteristics of the Algorithm.
C
C       SNSQE is a modification of the Powell hybrid method.  Two of
C       its main characteristics involve the choice of the correction as
C       a convex combination of the Newton and scaled gradient direc-
C       tions, and the updating of the Jacobian by the rank-1 method of
C       Broyden.  The choice of the correction guarantees (under reason-
C       able conditions) global convergence for starting points far from
C       the solution and a fast rate of convergence.  The Jacobian is
C       calculated at the starting point by either the user-supplied
C       subroutine or a forward-difference approximation, but it is not
C       recalculated until the rank-1 method fails to produce satis-
C       factory progress.
C
C       Timing.  The time required by SNSQE to solve a given problem
C         depends on N, the behavior of the functions, the accuracy
C         requested, and the starting point.  The number of arithmetic
C         operations needed by SNSQE is about 11.5*(N**2) to process
C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
C         to process each evaluation of the Jacobian (call to JAC,
C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
C         the timing of SNSQE will be strongly influenced by the time
C         spent in FCN and JAC.
C
C       Storage.  SNSQE requires (3*N**2 + 17*N)/2 single precision
C         storage locations, in addition to the storage required by the
C         program.  There are no internally declared storage arrays.
C
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), ..., X(9),
C       which solve the system of tridiagonal equations
C
C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
C                                   -X(8) + (3-2*X(9))*X(9) = -1
C
C       **********
C
C       PROGRAM TEST
C C
C C     Driver for SNSQE example.
C C
C       INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE
C       REAL TOL,FNORM
C       REAL X(9),FVEC(9),WA(180)
C       REAL ENORM,R1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 2
C       N = 9
C C
C C     The following starting values provide a rough solution.
C C
C       DO 10 J = 1, 9
C          X(J) = -1.E0
C    10    CONTINUE
C
C       LWA = 180
C       NPRINT = 0
C C
C C     Set TOL to the square root of the machine precision.
C C     Unless high precision solutions are required,
C C     this is the recommended setting.
C C
C       TOL = SQRT(R1MACH(4))
C C
C       CALL SNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
C       FNORM = ENORM(N,FVEC)
C       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
C       END
C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
C       INTEGER N,IFLAG
C       REAL X(N),FVEC(N)
C       INTEGER K
C       REAL ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C C
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
C       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
C
C       EXIT PARAMETER                         1
C
C       FINAL APPROXIMATE SOLUTION
C
C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
C
C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
C                 tions. In Numerical Methods for Nonlinear Algebraic
C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
C                 1988.
C***ROUTINES CALLED  SNSQ, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SNSQE
      INTEGER IOPT,N,NPRINT,INFO,LWA
      REAL TOL
      REAL X(*),FVEC(*),WA(LWA)
      EXTERNAL FCN, JAC
      INTEGER INDEX,J,LR,MAXFEV,ML,MODE,MU,NFEV,NJEV
      REAL EPSFCN,FACTOR,ONE,XTOL,ZERO
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0E2,1.0E0,0.0E0/
C***FIRST EXECUTABLE STATEMENT  SNSQE
      INFO = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0
     1    .OR. TOL .LT. ZERO .OR. LWA .LT. (3*N**2 +13*N)/2)
     2   GO TO 20
C
C     CALL SNSQ.
C
      MAXFEV = 100*(N + 1)
      IF (IOPT .EQ. 2) MAXFEV = 2*MAXFEV
      XTOL = TOL
      ML = N - 1
      MU = N - 1
      EPSFCN = ZERO
      MODE = 2
      DO 10 J = 1, N
         WA(J) = ONE
   10    CONTINUE
      LR = (N*(N + 1))/2
      INDEX=6*N+LR
      CALL SNSQ(FCN,JAC,IOPT,N,X,FVEC,WA(INDEX+1),N,XTOL,MAXFEV,ML,MU,
     1           EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
     2           WA(6*N+1),LR,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),
     3           WA(5*N+1))
      IF (INFO .EQ. 5) INFO = 4
   20 CONTINUE
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'SNSQE',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE SNSQE.
C
      END
