*DECK CPOFS
      SUBROUTINE CPOFS (A, LDA, N, V, ITASK, IND, WORK)
C***BEGIN PROLOGUE  CPOFS
C***PURPOSE  Solve a positive definite symmetric complex system of
C            linear equations.
C***LIBRARY   SLATEC
C***CATEGORY  D2D1B
C***TYPE      COMPLEX (SPOFS-S, DPOFS-D, CPOFS-C)
C***KEYWORDS  HERMITIAN, LINEAR EQUATIONS, POSITIVE DEFINITE, SYMMETRIC
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C    Subroutine CPOFS solves a  positive definite symmetric
C    NxN system of complex linear equations using LINPACK
C    subroutines CPOCO and CPOSL.  That is, if A is an NxN
C    complex positive definite symmetric matrix and if X and B
C    are complex N-vectors, then CPOFS solves the equation
C
C                          A*X=B.
C
C    Care should be taken not to use CPOFS with a non-Hermitian
C    matrix.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices R and R-TRANSPOSE.  These factors are used to
C    find the solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of a does not need to be performed again and
C    the option to only solve (ITASK .GT. 1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, and N must not have been altered by the user following
C    factorization (ITASK=1).  IND will not be changed by CPOFS
C    in this case.
C
C  Argument Description ***
C
C    A      COMPLEX(LDA,N)
C             on entry, the doubly subscripted array with dimension
C               (LDA,N) which contains the coefficient matrix.  Only
C               the upper triangle, including the diagonal, of the
C               coefficient matrix need be entered and will subse-
C               quently be referenced and changed by the routine.
C             on return, contains in its upper triangle an upper
C               triangular matrix R such that  A = (R-TRANSPOSE) * R .
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  N must be greater
C             than or equal to 1.  (terminal error message IND=-2)
C    V      COMPLEX(N)
C             on entry the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             if ITASK = 1, the matrix A is factored and then the
C               linear equation is solved.
C             if ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A.
C             if ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  see error message corresponding to IND below.
C    WORK   COMPLEX(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than 1.
C    IND=-3  terminal   ITASK is less than 1.
C    IND=-4  terminal   The matrix A is computationally singular or
C                         is not positive definite.  A solution
C                         has not been computed.
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the
C                         matrix A may be poorly scaled.
C
C               NOTE-  The above terminal(*fatal*) error messages are
C                      designed to be handled by XERMSG in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERMSG.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CPOCO, CPOSL, R1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800516  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to
C           IF-THEN-ELSE.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CPOFS
C
      INTEGER LDA,N,ITASK,IND,INFO
      COMPLEX A(LDA,*),V(*),WORK(*)
      REAL R1MACH
      REAL RCOND
      CHARACTER*8 XERN1, XERN2
C***FIRST EXECUTABLE STATEMENT  CPOFS
      IF (LDA.LT.N) THEN
         IND = -1
         WRITE (XERN1, '(I8)') LDA
         WRITE (XERN2, '(I8)') N
         CALL XERMSG ('SLATEC', 'CPOFS', 'LDA = ' // XERN1 //
     *      ' IS LESS THAN N = ' // XERN2, -1, 1)
         RETURN
      ENDIF
C
      IF (N.LE.0) THEN
         IND = -2
         WRITE (XERN1, '(I8)') N
         CALL XERMSG ('SLATEC', 'CPOFS', 'N = ' // XERN1 //
     *      ' IS LESS THAN 1', -2, 1)
         RETURN
      ENDIF
C
      IF (ITASK.LT.1) THEN
         IND = -3
         WRITE (XERN1, '(I8)') ITASK
         CALL XERMSG ('SLATEC', 'CPOFS', 'ITASK = ' // XERN1 //
     *      ' IS LESS THAN 1', -3, 1)
         RETURN
      ENDIF
C
      IF (ITASK.EQ.1) THEN
C
C        FACTOR MATRIX A INTO R
C
         CALL CPOCO(A,LDA,N,RCOND,WORK,INFO)
C
C        CHECK FOR POSITIVE DEFINITE MATRIX
C
         IF (INFO.NE.0) THEN
            IND = -4
            CALL XERMSG ('SLATEC', 'CPOFS',
     *         'SINGULAR OR NOT POSITIVE DEFINITE - NO SOLUTION', -4, 1)
            RETURN
         ENDIF
C
C        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
C        AND CHECK FOR IND GREATER THAN ZERO
C
         IND = -LOG10(R1MACH(4)/RCOND)
         IF (IND.LE.0) THEN
            IND = -10
            CALL XERMSG ('SLATEC', 'CPOFS',
     *         'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
         ENDIF
      ENDIF
C
C     SOLVE AFTER FACTORING
C
      CALL CPOSL(A,LDA,N,V)
      RETURN
      END
