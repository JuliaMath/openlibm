*DECK DPBSL
      SUBROUTINE DPBSL (ABD, LDA, N, M, B)
C***BEGIN PROLOGUE  DPBSL
C***PURPOSE  Solve a real symmetric positive definite band system
C            using the factors computed by DPBCO or DPBFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B2
C***TYPE      DOUBLE PRECISION (SPBSL-S, DPBSL-D, CPBSL-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX,
C             POSITIVE DEFINITE, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DPBSL solves the double precision symmetric positive definite
C     band system  A*X = B
C     using the factors computed by DPBCO or DPBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DPBCO or DPBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        M       INTEGER
C                the number of diagonals above the main diagonal.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal.  Technically this indicates
C        singularity, but it is usually caused by improper subroutine
C        arguments.  It will not occur if the subroutines are called
C        correctly, and  INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)
C           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DPBSL(ABD,LDA,N,C(1,J))
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPBSL
      INTEGER LDA,N,M
      DOUBLE PRECISION ABD(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,LA,LB,LM
C
C     SOLVE TRANS(R)*Y = B
C
C***FIRST EXECUTABLE STATEMENT  DPBSL
      DO 10 K = 1, N
         LM = MIN(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         T = DDOT(LM,ABD(LA,K),1,B(LB),1)
         B(K) = (B(K) - T)/ABD(M+1,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         LM = MIN(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         B(K) = B(K)/ABD(M+1,K)
         T = -B(K)
         CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
      RETURN
      END
