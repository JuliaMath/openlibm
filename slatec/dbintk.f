*DECK DBINTK
      SUBROUTINE DBINTK (X, Y, T, N, K, BCOEF, Q, WORK)
C***BEGIN PROLOGUE  DBINTK
C***PURPOSE  Compute the B-representation of a spline which interpolates
C            given data.
C***LIBRARY   SLATEC
C***CATEGORY  E1A
C***TYPE      DOUBLE PRECISION (BINTK-S, DBINTK-D)
C***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Abstract    **** a double precision routine ****
C
C         DBINTK is the SPLINT routine of the reference.
C
C         DBINTK produces the B-spline coefficients, BCOEF, of the
C         B-spline of order K with knots T(I), I=1,...,N+K, which
C         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
C         any of its derivatives can be evaluated by calls to DBVALU.
C
C         The I-th equation of the linear system A*BCOEF = B for the
C         coefficients of the interpolant enforces interpolation at
C         X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is
C         a band matrix with 2K-1 bands if A is invertible.  The matrix
C         A is generated row by row and stored, diagonal by diagonal,
C         in the rows of Q, with the main diagonal going into row K.
C         The banded system is then solved by a call to DBNFAC (which
C         constructs the triangular factorization for A and stores it
C         again in Q), followed by a call to DBNSLV (which then
C         obtains the solution BCOEF by substitution).  DBNFAC does no
C         pivoting, since the total positivity of the matrix A makes
C         this unnecessary.  The linear system to be solved is
C         (theoretically) invertible if and only if
C                 T(I) .LT. X(I) .LT. T(I+K),        for all I.
C         Equality is permitted on the left for I=1 and on the right
C         for I=N when K knots are used at X(1) or X(N).  Otherwise,
C         violation of this condition is certain to lead to an error.
C
C     Description of Arguments
C
C         Input       X,Y,T are double precision
C           X       - vector of length N containing data point abscissa
C                     in strictly increasing order.
C           Y       - corresponding vector of length N containing data
C                     point ordinates.
C           T       - knot vector of length N+K
C                     Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
C                     .GE. X(N), this leaves only N-K knots (not nec-
C                     essarily X(I) values) interior to (X(1),X(N))
C           N       - number of data points, N .GE. K
C           K       - order of the spline, K .GE. 1
C
C         Output      BCOEF,Q,WORK are double precision
C           BCOEF   - a vector of length N containing the B-spline
C                     coefficients
C           Q       - a work vector of length (2*K-1)*N, containing
C                     the triangular factorization of the coefficient
C                     matrix of the linear system being solved.  The
C                     coefficients for the interpolant of an
C                     additional data set (X(I),YY(I)), I=1,...,N
C                     with the same abscissa can be obtained by loading
C                     YY into BCOEF and then executing
C                         CALL DBNSLV (Q,2K-1,N,K-1,K-1,BCOEF)
C           WORK    - work vector of length 2*K
C
C     Error Conditions
C         Improper input is a fatal error
C         Singular system of equations is a fatal error
C
C***REFERENCES  D. E. Amos, Computation with splines and B-splines,
C                 Report SAND78-1968, Sandia Laboratories, March 1979.
C               Carl de Boor, Package for calculating with B-splines,
C                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
C                 pp. 441-472.
C               Carl de Boor, A Practical Guide to Splines, Applied
C                 Mathematics Series 27, Springer-Verlag, New York,
C                 1978.
C***ROUTINES CALLED  DBNFAC, DBNSLV, DBSPVN, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800901  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBINTK
C
      INTEGER IFLAG, IWORK, K, N, I, ILP1MX, J, JJ, KM1, KPKM2, LEFT,
     1 LENQ, NP1
      DOUBLE PRECISION BCOEF(*), Y(*), Q(*), T(*), X(*), XI, WORK(*)
C     DIMENSION Q(2*K-1,N), T(N+K)
C***FIRST EXECUTABLE STATEMENT  DBINTK
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      JJ = N - 1
      IF(JJ.EQ.0) GO TO 6
      DO 5 I=1,JJ
      IF(X(I).GE.X(I+1)) GO TO 110
    5 CONTINUE
    6 CONTINUE
      NP1 = N + 1
      KM1 = K - 1
      KPKM2 = 2*KM1
      LEFT = K
C                ZERO OUT ALL ENTRIES OF Q
      LENQ = N*(K+KM1)
      DO 10 I=1,LENQ
        Q(I) = 0.0D0
   10 CONTINUE
C
C  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN(I+K,NP1)
C        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
C                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
C        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
C        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
C        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
C        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
C        ARE RETURNED, IN  BCOEF (USED FOR TEMP. STORAGE HERE), BY THE
C        FOLLOWING
   30   CALL DBSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
C        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
C        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
C        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
C        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
C        DBNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
C        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
C        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
C        ENTRY
C            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
C                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
C        OF  Q .
        JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
        DO 40 J=1,K
          JJ = JJ + KPKM2
          Q(JJ) = BCOEF(J)
   40   CONTINUE
   50 CONTINUE
C
C     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL DBNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
C     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL DBNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
      RETURN
C
C
   80 CONTINUE
      CALL XERMSG ('SLATEC', 'DBINTK',
     +   'SOME ABSCISSA WAS NOT IN THE SUPPORT OF THE CORRESPONDING ' //
     +   'BASIS FUNCTION AND THE SYSTEM IS SINGULAR.', 2, 1)
      RETURN
   90 CONTINUE
      CALL XERMSG ('SLATEC', 'DBINTK',
     +   'THE SYSTEM OF SOLVER DETECTS A SINGULAR SYSTEM ALTHOUGH ' //
     +   'THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATISFIED.',
     +   8, 1)
      RETURN
  100 CONTINUE
      CALL XERMSG ('SLATEC', 'DBINTK', 'K DOES NOT SATISFY K.GE.1', 2,
     +   1)
      RETURN
  105 CONTINUE
      CALL XERMSG ('SLATEC', 'DBINTK', 'N DOES NOT SATISFY N.GE.K', 2,
     +   1)
      RETURN
  110 CONTINUE
      CALL XERMSG ('SLATEC', 'DBINTK',
     +   'X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR SOME I', 2, 1)
      RETURN
      END
