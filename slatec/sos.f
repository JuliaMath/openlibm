*DECK SOS
      SUBROUTINE SOS (FNC, NEQ, X, RTOLX, ATOLX, TOLF, IFLAG, RW, LRW,
     +   IW, LIW)
C***BEGIN PROLOGUE  SOS
C***PURPOSE  Solve a square system of nonlinear equations.
C***LIBRARY   SLATEC
C***CATEGORY  F2A
C***TYPE      SINGLE PRECISION (SOS-S, DSOS-D)
C***KEYWORDS  BROWN'S METHOD, NEWTON'S METHOD, NONLINEAR EQUATIONS,
C             ROOTS, SOLUTIONS
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C     SOS solves a system of NEQ simultaneous nonlinear equations in
C     NEQ unknowns.  That is, it solves the problem   F(X)=0
C     where X is a vector with components  X(1),...,X(NEQ)  and  F
C     is a vector of nonlinear functions.  Each equation is of the form
C
C               F (X(1),...,X(NEQ))=0     for K=1,...,NEQ.
C                K
C
C     The algorithm is based on an iterative method which is a
C     variation of Newton's method using Gaussian elimination
C     in a manner similar to the Gauss-Seidel process.  Convergence
C     is roughly quadratic.  All partial derivatives required by
C     the algorithm are approximated by first difference quotients.
C     The convergence behavior of this code is affected by the
C     ordering of the equations, and it is advantageous to place linear
C     and mildly nonlinear equations first in the ordering.
C
C     Actually, SOS is merely an interfacing routine for
C     calling subroutine SOSEQS which embodies the solution
C     algorithm.  The purpose of this is to add greater
C     flexibility and ease of use for the prospective user.
C
C     SOSEQS calls the accompanying routine SOSSOL, which solves special
C     triangular linear systems by back-substitution.
C
C     The user must supply a function subprogram which evaluates the
C     K-th equation only (K specified by SOSEQS) for each call
C     to the subprogram.
C
C     SOS represents an implementation of the mathematical algorithm
C     described in the references below.  It is a modification of the
C     code SOSNLE written by H. A. Watts in 1973.
C
C **********************************************************************
C   -Input-
C
C     FNC -Name of the function program which evaluates the equations.
C          This name must be in an EXTERNAL statement in the calling
C          program.  The user must supply FNC in the form FNC(X,K),
C          where X is the solution vector (which must be dimensioned
C          in FNC) and FNC returns the value of the K-th function.
C
C     NEQ -Number of equations to be solved.
C
C     X   -Solution vector.  Initial guesses must be supplied.
C
C     RTOLX -Relative error tolerance used in the convergence criteria.
C          Each solution component X(I) is checked by an accuracy test
C          of the form   ABS(X(I)-XOLD(I)) .LE. RTOLX*ABS(X(I))+ATOLX,
C          where XOLD(I) represents the previous iteration value.
C          RTOLX must be non-negative.
C
C     ATOLX -Absolute error tolerance used in the convergence criteria.
C          ATOLX must be non-negative.  If the user suspects some
C          solution component may be zero, he should set ATOLX to an
C          appropriate (depends on the scale of the remaining variables)
C          positive value for better efficiency.
C
C     TOLF -Residual error tolerance used in the convergence criteria.
C          Convergence will be indicated if all residuals (values of the
C          functions or equations) are not bigger than TOLF in
C          magnitude.  Note that extreme care must be given in assigning
C          an appropriate value for TOLF because this convergence test
C          is dependent on the scaling of the equations.  An
C          inappropriate value can cause premature termination of the
C          iteration process.
C
C     IFLAG -Optional input indicator.  You must set  IFLAG=-1  if you
C          want to use any of the optional input items listed below.
C          Otherwise set it to zero.
C
C     RW  -A REAL work array which is split apart by SOS and used
C          internally by SOSEQS.
C
C     LRW -Dimension of the RW array.  LRW must be at least
C                    1 + 6*NEQ + NEQ*(NEQ+1)/2
C
C     IW  -An INTEGER work array which is split apart by SOS and used
C          internally by SOSEQS.
C
C     LIW -Dimension of the IW array. LIW must be at least  3 + NEQ.
C
C   -Optional Input-
C
C     IW(1) -Internal printing parameter.  You must set  IW(1)=-1  if
C          you want the intermediate solution iterates to be printed.
C
C     IW(2) -Iteration limit.  The maximum number of allowable
C          iterations can be specified, if desired.  To override the
C          default value of 50, set IW(2) to the number wanted.
C
C     Remember, if you tell the code that you are using one of the
C               options (by setting IFLAG=-1), you must supply values
C               for both IW(1) and IW(2).
C
C **********************************************************************
C   -Output-
C
C     X   -Solution vector.
C
C     IFLAG -Status indicator
C
C                         *** Convergence to a Solution ***
C
C          1 Means satisfactory convergence to a solution was achieved.
C            Each solution component X(I) satisfies the error tolerance
C            test   ABS(X(I)-XOLD(I)) .LE. RTOLX*ABS(X(I))+ATOLX.
C
C          2 Means procedure converged to a solution such that all
C            residuals are at most TOLF in magnitude,
C            ABS(FNC(X,I)) .LE. TOLF.
C
C          3 Means that conditions for both IFLAG=1 and IFLAG=2 hold.
C
C          4 Means possible numerical convergence.  Behavior indicates
C            limiting precision calculations as a result of user asking
C            for too much accuracy or else convergence is very slow.
C            Residual norms and solution increment norms have
C            remained roughly constant over several consecutive
C            iterations.
C
C                         *** Task Interrupted ***
C
C          5 Means the allowable number of iterations has been met
C            without obtaining a solution to the specified accuracy.
C            Very slow convergence may be indicated.  Examine the
C            approximate solution returned and see if the error
C            tolerances seem appropriate.
C
C          6 Means the allowable number of iterations has been met and
C            the iterative process does not appear to be converging.
C            A local minimum may have been encountered or there may be
C            limiting precision difficulties.
C
C          7 Means that the iterative scheme appears to be diverging.
C            Residual norms and solution increment norms have
C            increased over several consecutive iterations.
C
C                         *** Task Cannot Be Continued ***
C
C          8 Means that a Jacobian-related matrix was singular.
C
C          9 Means improper input parameters.
C
C          *** IFLAG should be examined after each call to   ***
C          *** SOS with the appropriate action being taken.  ***
C
C
C     RW(1) -Contains a norm of the residual.
C
C     IW(3) -Contains the number of iterations used by the process.
C
C **********************************************************************
C***REFERENCES  K. M. Brown, Solution of simultaneous nonlinear
C                 equations, Algorithm 316, Communications of the
C                 A.C.M. 10, (1967), pp. 728-729.
C               K. M. Brown, A quadratically convergent Newton-like
C                 method based upon Gaussian elimination, SIAM Journal
C                 on Numerical Analysis 6, (1969), pp. 560-569.
C***ROUTINES CALLED  SOSEQS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900510  Convert XERRWV calls to XERMSG calls, changed Prologue
C           comments to agree with DSOS.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SOS
      DIMENSION X(*), RW(*), IW(*)
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
      EXTERNAL FNC
C***FIRST EXECUTABLE STATEMENT  SOS
      INPFLG = IFLAG
C
C     CHECK FOR VALID INPUT
C
      IF (NEQ .LE. 0) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'SOS', 'THE NUMBER OF EQUATIONS ' //
     *      'MUST BE A POSITIVE INTEGER.  YOU HAVE CALLED THE ' //
     *      'CODE WITH NEQ = ' // XERN1, 1, 1)
         IFLAG = 9
      ENDIF
C
      IF (RTOLX .LT. 0.0D0 .OR. ATOLX .LT. 0.0D0) THEN
         WRITE (XERN3, '(1PE15.6)') ATOLX
         WRITE (XERN4, '(1PE15.6)') RTOLX
         CALL XERMSG ('SLATEC', 'SOS', 'THE ERROR TOLERANCES FOR ' //
     *      'THE SOLUTION ITERATES CANNOT BE NEGATIVE. YOU HAVE ' //
     *      'CALLED THE CODE WITH  RTOLX = ' // XERN3 //
     *      ' AND ATOLX = ' // XERN4,2, 1)
            IFLAG = 9
      ENDIF
C
      IF (TOLF .LT. 0.0D0) THEN
         WRITE (XERN3, '(1PE15.6)') TOLF
         CALL XERMSG ('SLATEC', 'SOS', 'THE RESIDUAL ERROR ' //
     *      'TOLERANCE MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE ' //
     *      'CODE WITH TOLF = ' // XERN3, 3, 1)
            IFLAG = 9
      ENDIF
C
      IPRINT = 0
      MXIT = 50
      IF (INPFLG .EQ. (-1)) THEN
         IF (IW(1) .EQ. (-1)) IPRINT = -1
         MXIT = IW(2)
         IF (MXIT .LE. 0) THEN
            WRITE (XERN1, '(I8)') MXIT
            CALL XERMSG ('SLATEC', 'SOS', 'YOU HAVE TOLD THE CODE ' //
     *         'TO USE OPTIONAL IN PUT ITEMS BY SETTING  IFLAG=-1. ' //
     *         'HOWEVER YOU HAVE CALLED THE CODE WITH THE MAXIMUM ' //
     *         'ALLOWABLE NUMBER OF ITERATIONS SET TO  IW(2) = ' //
     *         XERN1, 4, 1)
            IFLAG = 9
         ENDIF
      ENDIF
C
      NC = (NEQ*(NEQ+1))/2
      IF (LRW .LT. 1 + 6*NEQ + NC) THEN
         WRITE (XERN1, '(I8)') LRW
         CALL XERMSG ('SLATEC', 'SOS', 'DIMENSION OF THE RW ARRAY ' //
     *      'MUST BE AT LEAST 1 + 6*NEQ + NEQ*(NEQ+1)/2 .  YOU HAVE ' //
     *      'CALLED THE CODE WITH LRW = ' // XERN1, 5, 1)
         IFLAG = 9
      ENDIF
C
      IF (LIW .LT. 3 + NEQ) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'SOS', 'DIMENSION OF THE IW ARRAY ' //
     *      'MUST BE AT LEAST   3 + NEQ.  YOU HAVE CALLED THE CODE ' //
     *      'WITH  LIW = ' // XERN1, 6, 1)
         IFLAG = 9
      ENDIF
C
      IF (IFLAG .NE. 9) THEN
         NCJS = 6
         NSRRC = 4
         NSRI = 5
C
         K1 = NC + 2
         K2 = K1 + NEQ
         K3 = K2 + NEQ
         K4 = K3 + NEQ
         K5 = K4 + NEQ
         K6 = K5 + NEQ
C
         CALL SOSEQS(FNC, NEQ, X, RTOLX, ATOLX, TOLF, IFLAG, MXIT, NCJS,
     1               NSRRC, NSRI, IPRINT, RW(1), RW(2), NC, RW(K1),
     2               RW(K2), RW(K3), RW(K4), RW(K5), RW(K6), IW(4))
C
         IW(3) = MXIT
      ENDIF
      RETURN
      END
