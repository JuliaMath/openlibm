*DECK SNLS1E
      SUBROUTINE SNLS1E (FCN, IOPT, M, N, X, FVEC, TOL, NPRINT, INFO,
     +   IW, WA, LWA)
C***BEGIN PROLOGUE  SNLS1E
C***PURPOSE  An easy-to-use code which minimizes the sum of the squares
C            of M nonlinear functions in N variables by a modification
C            of the Levenberg-Marquardt algorithm.
C***LIBRARY   SLATEC
C***CATEGORY  K1B1A1, K1B1A2
C***TYPE      SINGLE PRECISION (SNLS1E-S, DNLS1E-D)
C***KEYWORDS  EASY-TO-USE, LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
C             NONLINEAR LEAST SQUARES
C***AUTHOR  Hiebert, K. L., (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of SNLS1E is to minimize the sum of the squares of M
C       nonlinear functions in N variables by a modification of the
C       Levenberg-Marquardt algorithm.  This is done by using the more
C       general least-squares solver SNLS1.  The user must provide a
C       subroutine which calculates the functions.  The user has the
C       option of how the Jacobian will be supplied.  The user can
C       supply the full Jacobian, or the rows of the Jacobian (to avoid
C       storing the full Jacobian), or let the code approximate the
C       Jacobian by forward-differencing.  This code is the combination
C       of the MINPACK codes (Argonne) LMDER1, LMDIF1, and LMSTR1.
C
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE SNLS1E(FCN,IOPT,M,N,X,FVEC,TOL,NPRINT,
C      *                  INFO,IW,WA,LWA)
C       INTEGER IOPT,M,N,NPRINT,INFO,LWA
C       INTEGER IW(N)
C       REAL TOL
C       REAL X(N),FVEC(M),WA(LWA)
C       EXTERNAL FCN
C
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to SNLS1E and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from SNLS1E.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  If the user wants to supply the Jacobian
C         (IOPT=2 or 3), then FCN must be written to calculate the
C         Jacobian, as well as the functions.  See the explanation
C         of the IOPT argument below.
C         If the user wants the iterates printed (NPRINT positive), then
C         FCN must do the printing.  See the explanation of NPRINT
C         below.  FCN must be declared in an EXTERNAL statement in the
C         calling program and should be written as follows.
C
C
C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C         INTEGER IFLAG,LDFJAC,M,N
C         REAL X(N),FVEC(M)
C         ----------
C         FJAC and LDFJAC may be ignored     , if IOPT=1.
C         REAL FJAC(LDFJAC,N)                , if IOPT=2.
C         REAL FJAC(N)                       , if IOPT=3.
C         ----------
C           If IFLAG=0, the values in X and FVEC are available
C           for printing.  See the explanation of NPRINT below.
C           IFLAG will never be zero unless NPRINT is positive.
C           The values of X and FVEC must not be changed.
C         RETURN
C         ----------
C           If IFLAG=1, calculate the functions at X and return
C           this vector in FVEC.
C         RETURN
C         ----------
C           If IFLAG=2, calculate the full Jacobian at X and return
C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
C           IOPT=2.  FVEC contains the function values at X and must
C           not be altered.  FJAC(I,J) must be set to the derivative
C           of FVEC(I) with respect to X(J).
C         RETURN
C         ----------
C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
C           and return this vector in FJAC.  Note that IFLAG will
C           never be 3 unless IOPT=3.  FVEC contains the function
C           values at X and must not be altered.  FJAC(J) must be
C           set to the derivative of FVEC(LDFJAC) with respect to X(J).
C         RETURN
C         ----------
C         END
C
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of SNLS1E.  In this case,
C         set IFLAG to a negative integer.
C
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=2 or 3, then the user must supply the
C         Jacobian, as well as the function values, through the
C         subroutine FCN.  If IOPT=2, the user supplies the full
C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
C         one row of the Jacobian with each call.  (In this manner,
C         storage can be saved because the full Jacobian is not stored.)
C         If IOPT=1, the code will approximate the Jacobian by forward
C         differencing.
C
C       M is a positive integer input variable set to the number of
C         functions.
C
C       N is a positive integer input variable set to the number of
C         variables.  N must not exceed M.
C
C       X is an array of length N.  On input, X must contain an initial
C         estimate of the solution vector.  On output, X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length M which contains the functions
C         evaluated at the output X.
C
C       TOL is a non-negative input variable.  Termination occurs when
C         the algorithm estimates either that the relative error in the
C         sum of squares is at most TOL or that the relative error
C         between X and the solution is at most TOL.  Section 4 contains
C         more details about TOL.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. Appropriate
C         print statements must be added to FCN (see example) and
C         FVEC should not be altered.  If NPRINT is not positive, no
C         special calls of FCN with IFLAG = 0 are made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  improper input parameters.
C
C         INFO = 1  algorithm estimates that the relative error in the
C                   sum of squares is at most TOL.
C
C         INFO = 2  algorithm estimates that the relative error between
C                   X and the solution is at most TOL.
C
C         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
C
C         INFO = 4  FVEC is orthogonal to the columns of the Jacobian to
C                   machine precision.
C
C         INFO = 5  number of calls to FCN has reached 100*(N+1)
C                   for IOPT=2 or 3 or 200*(N+1) for IOPT=1.
C
C         INFO = 6  TOL is too small.  No further reduction in the sum
C                   of squares is possible.
C
C         INFO = 7  TOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       IW is an INTEGER work array of length N.
C
C       WA is a work array of length LWA.
C
C       LWA is a positive integer input variable not less than
C         N*(M+5)+M for IOPT=1 and 2 or N*(N+5)+M for IOPT=3.
C
C
C 4. Successful Completion.
C
C       The accuracy of SNLS1E is controlled by the convergence parame-
C       ter TOL.  This parameter is used in tests which make three types
C       of comparisons between the approximation X and a solution XSOL.
C       SNLS1E terminates when any of the tests is satisfied.  If TOL is
C       less than the machine precision (as defined by the function
C       R1MACH(4)), then SNLS1E only attempts to satisfy the test
C       defined by the machine precision.  Further progress is not usu-
C       ally possible.  Unless high precision solutions are required,
C       the recommended value for TOL is the square root of the machine
C       precision.
C
C       The tests assume that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently.  If these conditions
C       are not satisfied, then SNLS1E may incorrectly indicate conver-
C       gence.  If the Jacobian is coded correctly or IOPT=1,
C       then the validity of the answer can be checked, for example, by
C       rerunning SNLS1E with tighter tolerances.
C
C       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
C         of a vector Z, then this test attempts to guarantee that
C
C               ENORM(FVEC) .LE. (1+TOL)*ENORM(FVECS),
C
C         where FVECS denotes the functions evaluated at XSOL.  If this
C         condition is satisfied with TOL = 10**(-K), then the final
C         residual norm ENORM(FVEC) has K significant decimal digits and
C         INFO is set to 1 (or to 3 if the second test is also satis-
C         fied).
C
C       Second Convergence Test.  If D is a diagonal matrix (implicitly
C         generated by SNLS1E) whose entries contain scale factors for
C         the variables, then this test attempts to guarantee that
C
C               ENORM(D*(X-XSOL)) .LE.  TOL*ENORM(D*XSOL).
C
C         If this condition is satisfied with TOL = 10**(-K), then the
C         larger components of D*X have K significant decimal digits and
C         INFO is set to 2 (or to 3 if the first test is also satis-
C         fied).  There is a danger that the smaller components of D*X
C         may have large relative errors, but the choice of D is such
C         that the accuracy of the components of X is usually related to
C         their sensitivity.
C
C       Third Convergence Test.  This test is satisfied when FVEC is
C         orthogonal to the columns of the Jacobian to machine preci-
C         sion.  There is no clear relationship between this test and
C         the accuracy of SNLS1E, and furthermore, the test is equally
C         well satisfied at other critical points, namely maximizers and
C         saddle points.  Therefore, termination caused by this test
C         (INFO = 4) should be examined carefully.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of SNLS1E can be due to improper input
C       parameters, arithmetic interrupts, or an excessive number of
C       function evaluations.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1
C         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or TOL .LT. 0.E0,
C         or for IOPT=1 or 2 LWA .LT. N*(M+5)+M, or for IOPT=3
C         LWA .LT. N*(N+5)+M.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by SNLS1E.  In this
C         case, it may be possible to remedy the situation by not evalu-
C         ating the functions here, but instead setting the components
C         of FVEC to numbers that exceed those in the initial FVEC.
C
C       Excessive Number of Function Evaluations.  If the number of
C         calls to FCN reaches 100*(N+1) for IOPT=2 or 3 or 200*(N+1)
C         for IOPT=1, then this indicates that the routine is converging
C         very slowly as measured by the progress of FVEC, and INFO is
C         set to 5.  In this case, it may be helpful to restart SNLS1E,
C         thereby forcing it to disregard old (and possibly harmful)
C         information.
C
C
C 6. Characteristics of the Algorithm.
C
C       SNLS1E is a modification of the Levenberg-Marquardt algorithm.
C       Two of its main characteristics involve the proper use of
C       implicitly scaled variables and an optimal choice for the cor-
C       rection.  The use of implicitly scaled variables achieves scale
C       invariance of SNLS1E and limits the size of the correction in
C       any direction where the functions are changing rapidly.  The
C       optimal choice of the correction guarantees (under reasonable
C       conditions) global convergence from starting points far from the
C       solution and a fast rate of convergence for problems with small
C       residuals.
C
C       Timing.  The time required by SNLS1E to solve a given problem
C         depends on M and N, the behavior of the functions, the accu-
C         racy requested, and the starting point.  The number of arith-
C         metic operations needed by SNLS1E is about N**3 to process
C         each evaluation of the functions (call to FCN) and to process
C         each evaluation of the Jacobian SNLS1E takes M*N**2 for IOPT=2
C         (one call to JAC), M*N**2 for IOPT=1 (N calls to FCN) and
C         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
C         can be evaluated quickly, the timing of SNLS1E will be
C         strongly influenced by the time spent in FCN.
C
C       Storage.  SNLS1E requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
C         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
C         locations and N integer storage locations, in addition to
C         the storage required by the program.  There are no internally
C         declared storage arrays.
C
C *Long Description:
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), and X(3)
C       which provide the best fit (in the least squares sense) of
C
C             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
C
C       to the data
C
C             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
C                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
C
C       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
C       I-th component of FVEC is thus defined by
C
C             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
C
C       **********
C
C       PROGRAM TEST
C C
C C     Driver for SNLS1E example.
C C
C       INTEGER I,IOPT,M,N,NPRINT,JNFO,LWA,NWRITE
C       INTEGER IW(3)
C       REAL TOL,FNORM
C       REAL X(3),FVEC(15),WA(75)
C       REAL ENORM,R1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 1
C       M = 15
C       N = 3
C C
C C     The following starting values provide a rough fit.
C C
C       X(1) = 1.E0
C       X(2) = 1.E0
C       X(3) = 1.E0
C C
C       LWA = 75
C       NPRINT = 0
C C
C C     Set TOL to the square root of the machine precision.
C C     Unless high precision solutions are required,
C C     this is the recommended setting.
C C
C       TOL = SQRT(R1MACH(4))
C C
C       CALL SNLS1E(FCN,IOPT,M,N,X,FVEC,TOL,NPRINT,
C      *            INFO,IW,WA,LWA)
C       FNORM = ENORM(M,FVEC)
C       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
C       END
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
C C     This is the form of the FCN routine if IOPT=1,
C C     that is, if the user does not calculate the Jacobian.
C       INTEGER M,N,IFLAG
C       REAL X(N),FVEC(M)
C       INTEGER I
C       REAL TMP1,TMP2,TMP3,TMP4
C       REAL Y(15)
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C       END
C
C
C       Results obtained with different compilers or machines
C       may be slightly different.
C
C       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
C
C       EXIT PARAMETER                         1
C
C       FINAL APPROXIMATE SOLUTION
C
C        0.8241058E-01  0.1133037E+01  0.2343695E+01
C
C
C       For IOPT=2, FCN would be modified as follows to also
C       calculate the full Jacobian when IFLAG=2.
C
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C C
C C     This is the form of the FCN routine if IOPT=2,
C C     that is, if the user calculates the full Jacobian.
C C
C       INTEGER LDFJAC,M,N,IFLAG
C       REAL X(N),FVEC(M)
C       REAL FJAC(LDFJAC,N)
C       INTEGER I
C       REAL TMP1,TMP2,TMP3,TMP4
C       REAL Y(15)
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       IF(IFLAG.NE.1) GO TO 20
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C C
C C     Below, calculate the full Jacobian.
C C
C    20    CONTINUE
C C
C       DO 30 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
C          FJAC(I,1) = -1.E0
C          FJAC(I,2) = TMP1*TMP2/TMP4
C          FJAC(I,3) = TMP1*TMP3/TMP4
C    30    CONTINUE
C       RETURN
C       END
C
C
C       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
C         LDFJAC would be set to 3, and FCN would be written as
C         follows to calculate a row of the Jacobian when IFLAG=3.
C
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C C     This is the form of the FCN routine if IOPT=3,
C C     that is, if the user calculates the Jacobian row by row.
C       INTEGER M,N,IFLAG
C       REAL X(N),FVEC(M)
C       REAL FJAC(N)
C       INTEGER I
C       REAL TMP1,TMP2,TMP3,TMP4
C       REAL Y(15)
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       IF( IFLAG.NE.1) GO TO 20
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C C
C C     Below, calculate the LDFJAC-th row of the Jacobian.
C C
C    20 CONTINUE
C
C       I = LDFJAC
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
C          FJAC(1) = -1.E0
C          FJAC(2) = TMP1*TMP2/TMP4
C          FJAC(3) = TMP1*TMP3/TMP4
C       RETURN
C       END
C
C***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
C                 implementation and theory.  In Numerical Analysis
C                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
C                 Watson, Editor), Lecture Notes in Mathematics 630,
C                 Springer-Verlag, 1978.
C***ROUTINES CALLED  SNLS1, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890206  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SNLS1E
      INTEGER M,N,NPRINT,INFO,LWA,IOPT
      INTEGER INDEX,IW(*)
      REAL TOL
      REAL X(*),FVEC(*),WA(*)
      EXTERNAL FCN
      INTEGER MAXFEV,MODE,NFEV,NJEV
      REAL FACTOR,FTOL,GTOL,XTOL,ZERO,EPSFCN
      SAVE FACTOR, ZERO
      DATA FACTOR,ZERO /1.0E2,0.0E0/
C***FIRST EXECUTABLE STATEMENT  SNLS1E
      INFO = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3 .OR.
     1    N .LE. 0 .OR. M .LT. N .OR. TOL .LT. ZERO
     2    .OR. LWA .LT. N*(N+5) + M) GO TO 10
      IF (IOPT .LT. 3 .AND. LWA .LT. N*(M+5) + M) GO TO 10
C
C     CALL SNLS1.
C
      MAXFEV = 100*(N + 1)
      IF (IOPT .EQ. 1) MAXFEV = 2*MAXFEV
      FTOL = TOL
      XTOL = TOL
      GTOL = ZERO
      EPSFCN = ZERO
      MODE = 1
      INDEX = 5*N+M
      CALL SNLS1(FCN,IOPT,M,N,X,FVEC,WA(INDEX+1),M,FTOL,XTOL,GTOL,
     1           MAXFEV,EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
     2           IW,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))
      IF (INFO .EQ. 8) INFO = 4
   10 CONTINUE
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'SNLS1E',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE SNLS1E.
C
      END
