*DECK CGESL
      SUBROUTINE CGESL (A, LDA, N, IPVT, B, JOB)
C***BEGIN PROLOGUE  CGESL
C***PURPOSE  Solve the complex system A*X=B or CTRANS(A)*X=B using the
C            factors computed by CGECO or CGEFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1
C***TYPE      COMPLEX (SGESL-S, DGESL-D, CGESL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGESL solves the complex system
C     A * X = B  or  CTRANS(A) * X = B
C     using the factors computed by CGECO or CGEFA.
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the output from CGECO or CGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from CGECO or CGEFA.
C
C        B       COMPLEX(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  CTRANS(A)*X = B  where
C                            CTRANS(A)  is the conjugate transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if CGECO has set RCOND .GT. 0.0
C        or CGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL CGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL CGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CDOTC
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGESL
      INTEGER LDA,N,IPVT(*),JOB
      COMPLEX A(LDA,*),B(*)
C
      COMPLEX CDOTC,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  CGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL CAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL CAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
C        FIRST SOLVE  CTRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = CDOTC(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/CONJG(A(K,K))
   60    CONTINUE
C
C        NOW SOLVE CTRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + CDOTC(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
