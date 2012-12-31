*DECK CSIDI
      SUBROUTINE CSIDI (A, LDA, N, KPVT, DET, WORK, JOB)
C***BEGIN PROLOGUE  CSIDI
C***PURPOSE  Compute the determinant and inverse of a complex symmetric
C            matrix using the factors from CSIFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1, D3C1
C***TYPE      COMPLEX (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
C             SYMMETRIC
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     CSIDI computes the determinant and inverse
C     of a complex symmetric matrix using the factors from CSIFA.
C
C     On Entry
C
C        A       COMPLEX(LDA,N)
C                the output from CSIFA.
C
C        LDA     INTEGER
C                the leading dimension of the array A .
C
C        N       INTEGER
C                the order of the matrix A .
C
C        KVPT    INTEGER(N)
C                the pivot vector from CSIFA.
C
C        WORK    COMPLEX(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                JOB has the decimal expansion  AB  where
C                   If  B .NE. 0, the inverse is computed,
C                   If  A .NE. 0, the determinant is computed,
C
C                For example, JOB = 11  gives both.
C
C     On Return
C
C        Variables not requested by JOB are not used.
C
C        A      contains the upper triangle of the inverse of
C               the original matrix.  The strict lower triangle
C               is never referenced.
C
C        DET    COMPLEX(2)
C               determinant of original matrix.
C               Determinant = DET(1) * 10.0**DET(2)
C               with 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               or DET(1) = 0.0.
C
C     Error Condition
C
C        A division by zero may occur if the inverse is requested
C        and  CSICO  has set RCOND .EQ. 0.0
C        or  CSIFA  has set  INFO .NE. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CCOPY, CDOTU, CSWAP
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891107  Corrected category and modified routine equivalence
C           list.  (WRB)
C   891107  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSIDI
      INTEGER LDA,N,JOB
      COMPLEX A(LDA,*),DET(2),WORK(*)
      INTEGER KPVT(*)
C
      COMPLEX AK,AKP1,AKKP1,CDOTU,D,T,TEMP
      REAL TEN
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NOINV,NODET
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C***FIRST EXECUTABLE STATEMENT  CSIDI
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
C
      IF (NODET) GO TO 100
         DET(1) = (1.0E0,0.0E0)
         DET(2) = (0.0E0,0.0E0)
         TEN = 10.0E0
         T = (0.0E0,0.0E0)
         DO 90 K = 1, N
            D = A(K,K)
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 30
C
C              2 BY 2 BLOCK
C              USE DET (D  T)  =  (D/T * C - T) * T
C                      (T  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (CABS1(T) .NE. 0.0E0) GO TO 10
                  T = A(K,K+1)
                  D = (D/T)*A(K+1,K+1) - T
               GO TO 20
   10          CONTINUE
                  D = T
                  T = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
C
            DET(1) = D*DET(1)
            IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 80
   40          IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 50
                  DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                  DET(2) = DET(2) - (1.0E0,0.0E0)
               GO TO 40
   50          CONTINUE
   60          IF (CABS1(DET(1)) .LT. TEN) GO TO 70
                  DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                  DET(2) = DET(2) + (1.0E0,0.0E0)
               GO TO 60
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 230
         K = 1
  110    IF (K .GT. N) GO TO 220
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 140
C
C              1 BY 1
C
               A(K,K) = (1.0E0,0.0E0)/A(K,K)
               IF (KM1 .LT. 1) GO TO 130
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 120 J = 1, KM1
                     A(J,K) = CDOTU(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  120             CONTINUE
                  A(K,K) = A(K,K) + CDOTU(KM1,WORK,1,A(1,K),1)
  130          CONTINUE
               KSTEP = 1
            GO TO 180
  140       CONTINUE
C
C              2 BY 2
C
               T = A(K,K+1)
               AK = A(K,K)/T
               AKP1 = A(K+1,K+1)/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - (1.0E0,0.0E0))
               A(K,K) = AKP1/D
               A(K+1,K+1) = AK/D
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 170
                  CALL CCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 150 J = 1, KM1
                     A(J,K+1) = CDOTU(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  150             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1)
     1                         + CDOTU(KM1,WORK,1,A(1,K+1),1)
                  A(K,K+1) = A(K,K+1) + CDOTU(KM1,A(1,K),1,A(1,K+1),1)
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = CDOTU(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K) + CDOTU(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
               KSTEP = 2
  180       CONTINUE
C
C           SWAP
C
            KS = ABS(KPVT(K))
            IF (KS .EQ. K) GO TO 210
               CALL CSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 190 JB = KS, K
                  J = K + KS - JB
                  TEMP = A(J,K)
                  A(J,K) = A(KS,J)
                  A(KS,J) = TEMP
  190          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 200
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  200          CONTINUE
  210       CONTINUE
            K = K + KSTEP
         GO TO 110
  220    CONTINUE
  230 CONTINUE
      RETURN
      END
