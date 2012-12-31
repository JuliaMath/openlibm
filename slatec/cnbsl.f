*DECK CNBSL
      SUBROUTINE CNBSL (ABE, LDA, N, ML, MU, IPVT, B, JOB)
C***BEGIN PROLOGUE  CNBSL
C***PURPOSE  Solve a complex band system using the factors computed by
C            CNBCO or CNBFA.
C***LIBRARY   SLATEC
C***CATEGORY  D2C2
C***TYPE      COMPLEX (SNBSL-S, DNBSL-D, CNBSL-C)
C***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC, SOLVE
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C     CNBSL solves the complex band system
C     A * X = B  or  CTRANS(A) * X = B
C     using the factors computed by CNBCO or CNBFA.
C
C     On Entry
C
C        ABE     COMPLEX(LDA, NC)
C                the output from CNBCO or CNBFA.
C                NC must be .GE. 2*ML+MU+1 .
C
C        LDA     INTEGER
C                the leading dimension of the array  ABE .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from CNBCO or CNBFA.
C
C        B       COMPLEX(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B .
C                = nonzero   to solve  CTRANS(A)*X = B , where
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
C        setting of LDA.  It will not occur if the subroutines are
C        called correctly and if CNBCO has set RCOND .GT. 0.0
C        or CNBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL CNBCO(ABE,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C             CALL CNBSL(ABE,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CDOTC
C***REVISION HISTORY  (YYMMDD)
C   800730  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CNBSL
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      COMPLEX ABE(LDA,*),B(*)
C
      COMPLEX CDOTC,T
      INTEGER K,KB,L,LB,LDB,LM,M,MLM,NM1
C***FIRST EXECUTABLE STATEMENT  CNBSL
      M=MU+ML+1
      NM1=N-1
      LDB=1-LDA
      IF(JOB.NE.0)GO TO 50
C
C       JOB = 0 , SOLVE  A * X = B
C       FIRST SOLVE L*Y = B
C
        IF(ML.EQ.0)GO TO 30
        IF(NM1.LT.1)GO TO 30
          DO 20 K=1,NM1
            LM=MIN(ML,N-K)
            L=IPVT(K)
            T=B(L)
            IF(L.EQ.K)GO TO 10
              B(L)=B(K)
              B(K)=T
   10       CONTINUE
            MLM=ML-(LM-1)
            CALL CAXPY(LM,T,ABE(K+LM,MLM),LDB,B(K+1),1)
   20     CONTINUE
   30   CONTINUE
C
C       NOW SOLVE  U*X = Y
C
        DO 40 KB=1,N
          K=N+1-KB
          B(K)=B(K)/ABE(K,ML+1)
          LM=MIN(K,M)-1
          LB=K-LM
          T=-B(K)
          CALL CAXPY(LM,T,ABE(K-1,ML+2),LDB,B(LB),1)
   40   CONTINUE
      GO TO 100
   50 CONTINUE
C
C       JOB = NONZERO, SOLVE CTRANS(A) * X = B
C       FIRST SOLVE  CTRANS(U)*Y = B
C
        DO 60 K = 1, N
          LM = MIN(K,M) - 1
          LB = K - LM
          T = CDOTC(LM,ABE(K-1,ML+2),LDB,B(LB),1)
          B(K) = (B(K) - T)/CONJG(ABE(K,ML+1))
   60   CONTINUE
C
C       NOW SOLVE CTRANS(L)*X = Y
C
        IF (ML .EQ. 0) GO TO 90
        IF (NM1 .LT. 1) GO TO 90
          DO 80 KB = 1, NM1
            K = N - KB
            LM = MIN(ML,N-K)
            MLM = ML - (LM - 1)
            B(K) = B(K) + CDOTC(LM,ABE(K+LM,MLM),LDB,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
              T = B(L)
              B(L) = B(K)
              B(K) = T
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
      RETURN
      END
