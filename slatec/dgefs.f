*DECK DGEFS
      SUBROUTINE DGEFS (A, LDA, N, V, ITASK, IND, WORK, IWORK)
C***BEGIN PROLOGUE  DGEFS
C***PURPOSE  Solve a general system of linear equations.
C***LIBRARY   SLATEC
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFS-S, DGEFS-D, CGEFS-C)
C***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
C             GENERAL SYSTEM OF LINEAR EQUATIONS
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C    Subroutine DGEFS solves a general NxN system of double
C    precision linear equations using LINPACK subroutines DGECO
C    and DGESL.  That is, if A is an NxN double precision matrix
C    and if X and B are double precision N-vectors, then DGEFS
C    solves the equation
C
C                          A*X=B.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices U and L using partial pivoting.  These
C    factors and the pivoting information are used to find the
C    solution vector X.  An approximate condition number is
C    calculated to provide a rough estimate of the number of
C    digits of accuracy in the computed solution.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option to only solve (ITASK.GT.1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, N and IWORK must not have been altered by the user follow-
C    ing factorization (ITASK=1).  IND will not be changed by DGEFS
C    in this case.
C
C  Argument Description ***
C
C    A      DOUBLE PRECISION(LDA,N)
C             on entry, the doubly subscripted array with dimension
C               (LDA,N) which contains the coefficient matrix.
C             on return, an upper triangular matrix U and the
C               multipliers necessary to construct a matrix L
C               so that A=L*U.
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  The first N elements of
C             the array A are the elements of the first column of
C             the matrix A.  N must be greater than or equal to 1.
C             (terminal error message IND=-2)
C    V      DOUBLE PRECISION(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             If ITASK=1, the matrix A is factored and then the
C               linear equation is solved.
C             If ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A and IWORK.
C             If ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT. 0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.
C             LT. 0  see error message corresponding to IND below.
C    WORK   DOUBLE PRECISION(N)
C             a singly subscripted array of dimension at least N.
C    IWORK  INTEGER(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than 1.
C    IND=-3  terminal   ITASK is less than 1.
C    IND=-4  terminal   The matrix A is computationally singular.
C                         A solution has not been computed.
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the matrix
C                         A may be poorly scaled.
C
C               Note-  The above terminal(*fatal*) error messages are
C                      designed to be handled by XERMSG in which
C                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
C                      for warning error messages from XERMSG.  Unless
C                      the user provides otherwise, an error message
C                      will be printed followed by an abort.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  D1MACH, DGECO, DGESL, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800326  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFS
C
      INTEGER LDA,N,ITASK,IND,IWORK(*)
      DOUBLE PRECISION A(LDA,*),V(*),WORK(*),D1MACH
      DOUBLE PRECISION RCOND
      CHARACTER*8 XERN1, XERN2
C***FIRST EXECUTABLE STATEMENT  DGEFS
      IF (LDA.LT.N) THEN
         IND = -1
         WRITE (XERN1, '(I8)') LDA
         WRITE (XERN2, '(I8)') N
         CALL XERMSG ('SLATEC', 'DGEFS', 'LDA = ' // XERN1 //
     *      ' IS LESS THAN N = ' // XERN2, -1, 1)
         RETURN
      ENDIF
C
      IF (N.LE.0) THEN
         IND = -2
         WRITE (XERN1, '(I8)') N
         CALL XERMSG ('SLATEC', 'DGEFS', 'N = ' // XERN1 //
     *      ' IS LESS THAN 1', -2, 1)
         RETURN
      ENDIF
C
      IF (ITASK.LT.1) THEN
         IND = -3
         WRITE (XERN1, '(I8)') ITASK
         CALL XERMSG ('SLATEC', 'DGEFS', 'ITASK = ' // XERN1 //
     *      ' IS LESS THAN 1', -3, 1)
         RETURN
      ENDIF
C
      IF (ITASK.EQ.1) THEN
C
C        FACTOR MATRIX A INTO LU
C
         CALL DGECO(A,LDA,N,IWORK,RCOND,WORK)
C
C        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
C
         IF (RCOND.EQ.0.0D0) THEN
            IND = -4
            CALL XERMSG ('SLATEC', 'DGEFS',
     *         'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
            RETURN
         ENDIF
C
C        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
C        AND CHECK FOR IND GREATER THAN ZERO
C
         IND = -LOG10(D1MACH(4)/RCOND)
         IF (IND.LE.0) THEN
            IND=-10
            CALL XERMSG ('SLATEC', 'DGEFS',
     *         'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
         ENDIF
      ENDIF
C
C     SOLVE AFTER FACTORING
C
      CALL DGESL(A,LDA,N,IWORK,V,0)
      RETURN
      END
