*DECK SSPSL
      SUBROUTINE SSPSL (AP, N, KPVT, B)
C***BEGIN PROLOGUE  SSPSL
C***PURPOSE  Solve a real symmetric system using the factors obtained
C            from SSPFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1A
C***TYPE      SINGLE PRECISION (SSPSL-S, DSPSL-D, CHPSL-C, CSPSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED, SOLVE, SYMMETRIC
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     SSISL solves the real symmetric system
C     A * X = B
C     using the factors computed by SSPFA.
C
C     On Entry
C
C        AP      REAL(N*(N+1)/2)
C                the output from SSPFA.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        KPVT    INTEGER(N)
C                the pivot vector from SSPFA.
C
C        B       REAL(N)
C                the right hand side vector.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero may occur if  SSPCO  has set RCOND .EQ. 0.0
C        or  SSPFA  has set INFO .NE. 0  .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL SSPFA(AP,N,KPVT,INFO)
C           IF (INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL SSPSL(AP,N,KPVT,C(1,J))
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SAXPY, SDOT
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
C***END PROLOGUE  SSPSL
      INTEGER N,KPVT(*)
      REAL AP(*),B(*)
C
      REAL AK,AKM1,BK,BKM1,SDOT,DENOM,TEMP
      INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
C***FIRST EXECUTABLE STATEMENT  SSPSL
      K = N
      IK = (N*(N - 1))/2
   10 IF (K .EQ. 0) GO TO 80
         KK = IK + K
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
               CALL SAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
C
C           APPLY D INVERSE.
C
            B(K) = B(K)/AP(KK)
            K = K - 1
            IK = IK - K
         GO TO 70
   40    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            IKM1 = IK - (K - 1)
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
               CALL SAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
               CALL SAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
C
C           APPLY D INVERSE.
C
            KM1K = IK + K - 1
            KK = IK + K
            AK = AP(KK)/AP(KM1K)
            KM1KM1 = IKM1 + K - 1
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = B(K)/AP(KM1K)
            BKM1 = B(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            B(K) = (AKM1*BK - BKM1)/DENOM
            B(K-1) = (AK*BKM1 - BK)/DENOM
            K = K - 2
            IK = IK - (K + 1) - K
   70    CONTINUE
      GO TO 10
   80 CONTINUE
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
C
      K = 1
      IK = 0
   90 IF (K .GT. N) GO TO 160
         IF (KPVT(K) .LT. 0) GO TO 120
C
C           1 X 1 PIVOT BLOCK.
C
            IF (K .EQ. 1) GO TO 110
C
C              APPLY THE TRANSFORMATION.
C
               B(K) = B(K) + SDOT(K-1,AP(IK+1),1,B(1),1)
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
            IK = IK + K
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
               B(K) = B(K) + SDOT(K-1,AP(IK+1),1,B(1),1)
               IKP1 = IK + K
               B(K+1) = B(K+1) + SDOT(K-1,AP(IKP1+1),1,B(1),1)
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
            IK = IK + K + K + 1
            K = K + 2
  150    CONTINUE
      GO TO 90
  160 CONTINUE
      RETURN
      END
