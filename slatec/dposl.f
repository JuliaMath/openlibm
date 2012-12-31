*DECK DPOSL
      SUBROUTINE DPOSL (A, LDA, N, B)
C***BEGIN PROLOGUE  DPOSL
C***PURPOSE  Solve the real symmetric positive definite linear system
C            using the factors computed by DPOCO or DPOFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1B
C***TYPE      DOUBLE PRECISION (SPOSL-S, DPOSL-D, CPOSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DPOSL solves the double precision symmetric positive definite
C     system A * X = B
C     using the factors computed by DPOCO or DPOFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DPOCO or DPOFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
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
C        correctly and  INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)
C           IF (RCOND is too small .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DPOSL(A,LDA,N,C(1,J))
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPOSL
      INTEGER LDA,N
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB
C
C     SOLVE TRANS(R)*Y = B
C
C***FIRST EXECUTABLE STATEMENT  DPOSL
      DO 10 K = 1, N
         T = DDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END
