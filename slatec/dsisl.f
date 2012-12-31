*DECK DSISL
      SUBROUTINE DSISL (A, LDA, N, KPVT, B)
C***BEGIN PROLOGUE  DSISL
C***PURPOSE  Solve a real symmetric system using the factors obtained
C            from SSIFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1A
C***TYPE      DOUBLE PRECISION (SSISL-S, DSISL-D, CHISL-C, CSISL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, SYMMETRIC
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     DSISL solves the double precision symmetric system
C     A * X = B
C     using the factors computed by DSIFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA,N)
C                the output from DSIFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        KPVT    INTEGER(N)
C                the pivot vector from DSIFA.
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
C        A division by zero may occur if  DSICO  has set RCOND .EQ. 0.0
C        or  DSIFA  has set INFO .NE. 0  .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DSIFA(A,LDA,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DSISL(A,LDA,N,KPVT,C(1,J))
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891107  Modified routine equivalence list.  (WRB)
C   891107  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSISL
      INTEGER LDA,N,KPVT(*)
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
      INTEGER K,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
C***FIRST EXECUTABLE STATEMENT  DSISL
      K = N
   10 IF (K .EQ. 0) GO TO 80
         IF (KPVT(K) .LT. 0) GO TO 40
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 30
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 20
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
   20          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/A(K,K)
            K = K - 1
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 2) GO TO 60
               KP = ABS(KPVT(K))
               IF (KP .EQ. K - 1) GO TO 50
C
C                 INTERCHANGE.
C
                  TEMP = B(K-1)
                  B(K-1) = B(KP)
                  B(KP) = TEMP
   50          CONTINUE
C
C              APPLY THE TRANSFORMATION.
C
               CALL DAXPY(K-2,B(K),A(1,K),1,B(1),1)
               CALL DAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            AK = A(K,K)/A(K-1,K)
            AKM1 = A(K-1,K-1)/A(K-1,K)
            BK = B(K)/A(K-1,K)
            BKM1 = B(K-1)/A(K-1,K)
            DENOM = AK*AKM1 - 1.0D0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               KP = KPVT(K)
               IF (KP .EQ. K) GO TO 100
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
            K = K + 1
         GO TO 150
  120    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 140
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + DDOT(K-1,A(1,K),1,B(1),1)
               B(K+1) = B(K+1) + DDOT(K-1,A(1,K+1),1,B(1),1)
               KP = ABS(KPVT(K))
               IF (KP .EQ. K) GO TO 130
C
C                 INTERCHANGE.
C
                  TEMP = B(K)
                  B(K) = B(KP)
                  B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
