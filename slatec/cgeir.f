*DECK CGEIR
      SUBROUTINE CGEIR (A, LDA, N, V, ITASK, IND, WORK, IWORK)
C***BEGIN PROLOGUE  CGEIR
C***PURPOSE  Solve a general system of linear equations.  Iterative
C            refinement is used to obtain an error estimate.
C***LIBRARY   SLATEC
C***CATEGORY  D2C1
C***TYPE      COMPLEX (SGEIR-S, CGEIR-C)
C***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
C             GENERAL SYSTEM OF LINEAR EQUATIONS
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C    Subroutine CGEIR solves a general NxN system of complex
C    linear equations using LINPACK subroutines CGEFA and CGESL.
C    One pass of iterative refinement is used only to obtain an
C    estimate of the accuracy.  That is, if A is an NxN complex
C    matrix and if X and B are complex N-vectors, then CGEIR solves
C    the equation
C
C                          A*X=B.
C
C    The matrix A is first factored into upper and lower tri-
C    angular matrices U and L using partial pivoting.  These
C    factors and the pivoting information are used to calculate
C    the solution, X.  Then the residual vector is found and
C    used to calculate an estimate of the relative error, IND.
C    IND estimates the accuracy of the solution only when the
C    input matrix and the right hand side are represented
C    exactly in the computer and does not take into
C    account any errors in the input data.
C
C    If the equation A*X=B is to be solved for more than one vector
C    B, the factoring of A does not need to be performed again and
C    the option to only solve (ITASK .GT. 1) will be faster for
C    the succeeding solutions.  In this case, the contents of A,
C    LDA, N, WORK, and IWORK must not have been altered by the
C    user following factorization (ITASK=1).  IND will not be
C    changed by CGEIR in this case.
C
C  Argument Description ***
C
C    A      COMPLEX(LDA,N)
C             the doubly subscripted array with dimension (LDA,N)
C             which contains the coefficient matrix.  A is not
C             altered by the routine.
C    LDA    INTEGER
C             the leading dimension of the array A.  LDA must be great-
C             er than or equal to N.  (Terminal error message IND=-1)
C    N      INTEGER
C             the order of the matrix A.  The first N elements of
C             the array A are the elements of the first column of
C             matrix A.  N must be greater than or equal to 1.
C             (Terminal error message IND=-2)
C    V      COMPLEX(N)
C             on entry, the singly subscripted array(vector) of di-
C               mension N which contains the right hand side B of a
C               system of simultaneous linear equations A*X=B.
C             on return, V contains the solution vector, X .
C    ITASK  INTEGER
C             if ITASK=1, the matrix A is factored and then the
C               linear equation is solved.
C             if ITASK .GT. 1, the equation is solved using the existing
C               factored matrix A (stored in work).
C             if ITASK .LT. 1, then terminal error message IND=-3 is
C               printed.
C    IND    INTEGER
C             GT.0  IND is a rough estimate of the number of digits
C                     of accuracy in the solution, X.  IND=75 means
C                     that the solution vector X is zero.
C             LT.0  see error message corresponding to IND below.
C    WORK   COMPLEX(N*(N+1))
C             a singly subscripted array of dimension at least N*(N+1).
C    IWORK  INTEGER(N)
C             a singly subscripted array of dimension at least N.
C
C  Error Messages Printed ***
C
C    IND=-1  terminal   N is greater than LDA.
C    IND=-2  terminal   N is less than one.
C    IND=-3  terminal   ITASK is less than one.
C    IND=-4  terminal   The matrix A is computationally singular.
C                         A solution has not been computed.
C    IND=-10 warning    The solution has no apparent significance.
C                         The solution may be inaccurate or the matrix
C                         A may be poorly scaled.
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
C***ROUTINES CALLED  CCOPY, CDCDOT, CGEFA, CGESL, R1MACH, SCASUM, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800502  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to
C           IF-THEN-ELSE.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGEIR
C
      INTEGER LDA,N,ITASK,IND,IWORK(*),INFO,J
      COMPLEX A(LDA,*),V(*),WORK(N,*),CDCDOT
      REAL SCASUM,XNORM,DNORM,R1MACH
      CHARACTER*8 XERN1, XERN2
C***FIRST EXECUTABLE STATEMENT  CGEIR
      IF (LDA.LT.N) THEN
         IND = -1
         WRITE (XERN1, '(I8)') LDA
         WRITE (XERN2, '(I8)') N
         CALL XERMSG ('SLATEC', 'CGEIR', 'LDA = ' // XERN1 //
     *      ' IS LESS THAN N = ' // XERN2, -1, 1)
         RETURN
      ENDIF
C
      IF (N.LE.0) THEN
         IND = -2
         WRITE (XERN1, '(I8)') N
         CALL XERMSG ('SLATEC', 'CGEIR', 'N = ' // XERN1 //
     *      ' IS LESS THAN 1', -2, 1)
         RETURN
      ENDIF
C
      IF (ITASK.LT.1) THEN
         IND = -3
         WRITE (XERN1, '(I8)') ITASK
         CALL XERMSG ('SLATEC', 'CGEIR', 'ITASK = ' // XERN1 //
     *      ' IS LESS THAN 1', -3, 1)
         RETURN
      ENDIF
C
      IF (ITASK.EQ.1) THEN
C        MOVE MATRIX A TO WORK
         DO 10 J=1,N
            CALL CCOPY(N,A(1,J),1,WORK(1,J),1)
   10    CONTINUE
C
C        FACTOR MATRIX A INTO LU
C
         CALL CGEFA(WORK,N,N,IWORK,INFO)
C
C        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
C
         IF (INFO.NE.0) THEN
            IND = -4
            CALL XERMSG ('SLATEC', 'CGEIR',
     *         'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
            RETURN
         ENDIF
      ENDIF
C
C     SOLVE WHEN FACTORING COMPLETE
C     MOVE VECTOR B TO WORK
C
      CALL CCOPY(N,V(1),1,WORK(1,N+1),1)
      CALL CGESL(WORK,N,N,IWORK,V,0)
C
C     FORM NORM OF X0
C
      XNORM = SCASUM(N,V(1),1)
      IF (XNORM.EQ.0.0) THEN
         IND = 75
         RETURN
      ENDIF
C
C     COMPUTE  RESIDUAL
C
      DO 40 J=1,N
         WORK(J,N+1) = CDCDOT(N,-WORK(J,N+1),A(J,1),LDA,V,1)
   40 CONTINUE
C
C     SOLVE A*DELTA=R
C
      CALL CGESL(WORK,N,N,IWORK,WORK(1,N+1),0)
C
C     FORM NORM OF DELTA
C
      DNORM = SCASUM(N,WORK(1,N+1),1)
C
C     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
C     AND CHECK FOR IND GREATER THAN ZERO
C
      IND = -LOG10(MAX(R1MACH(4),DNORM/XNORM))
      IF (IND.LE.0) THEN
         IND = -10
         CALL XERMSG ('SLATEC', 'CGEIR',
     *      'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
      ENDIF
      RETURN
      END
