*DECK CSPDI
      SUBROUTINE CSPDI (AP, N, KPVT, DET, WORK, JOB)
C***BEGIN PROLOGUE  CSPDI
C***PURPOSE  Compute the determinant and inverse of a complex symmetric
C            matrix stored in packed form using the factors from CSPFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1, D3C1
C***TYPE      COMPLEX (SSPDI-S, DSPDI-D, CHPDI-C, CSPDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
C             PACKED, SYMMETRIC
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     CSPDI computes the determinant and inverse
C     of a complex symmetric matrix using the factors from CSPFA,
C     where the matrix is stored in packed form.
C
C     On Entry
C
C        AP      COMPLEX (N*(N+1)/2)
C                the output from CSPFA.
C
C        N       INTEGER
C                the order of the matrix A .
C
C        KVPT    INTEGER(N)
C                the pivot vector from CSPFA.
C
C        WORK    COMPLEX(N)
C                work vector.  Contents ignored.
C
C        JOB     INTEGER
C                JOB has the decimal expansion  AB  where
C                   if  B .NE. 0, the inverse is computed,
C                   if  A .NE. 0, the determinant is computed.
C
C                For example, JOB = 11  gives both.
C
C     On Return
C
C        Variables not requested by JOB are not used.
C
C        AP     contains the upper triangle of the inverse of
C               the original matrix, stored in packed form.
C               The columns of the upper triangle are stored
C               sequentially in a one-dimensional array.
C
C        DET    COMPLEX(2)
C               determinant of original matrix.
C               Determinant = DET(1) * 10.0**DET(2)
C               with 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               or DET(1) = 0.0.
C
C     Error Condition
C
C        A division by zero will occur if the inverse is requested
C        and  CSPCO  has set RCOND .EQ. 0.0
C        or  CSPFA  has set  INFO .NE. 0 .
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
C***END PROLOGUE  CSPDI
      INTEGER N,JOB
      COMPLEX AP(*),WORK(*),DET(2)
      INTEGER KPVT(*)
C
      COMPLEX AK,AKKP1,AKP1,CDOTU,D,T,TEMP
      REAL TEN
      INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
      INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
      LOGICAL NOINV,NODET
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C
C***FIRST EXECUTABLE STATEMENT  CSPDI
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
C
      IF (NODET) GO TO 110
         DET(1) = (1.0E0,0.0E0)
         DET(2) = (0.0E0,0.0E0)
         TEN = 10.0E0
         T = (0.0E0,0.0E0)
         IK = 0
         DO 100 K = 1, N
            KK = IK + K
            D = AP(KK)
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
                  IKP1 = IK + K
                  KKP1 = IKP1 + K
                  T = AP(KKP1)
                  D = (D/T)*AP(KKP1+1) - T
               GO TO 20
   10          CONTINUE
                  D = T
                  T = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
C
            IF (NODET) GO TO 90
               DET(1) = D*DET(1)
               IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 80
   40             IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 50
                     DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                     DET(2) = DET(2) - (1.0E0,0.0E0)
                  GO TO 40
   50             CONTINUE
   60             IF (CABS1(DET(1)) .LT. TEN) GO TO 70
                     DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                     DET(2) = DET(2) + (1.0E0,0.0E0)
                  GO TO 60
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
            IK = IK + K
  100    CONTINUE
  110 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 240
         K = 1
         IK = 0
  120    IF (K .GT. N) GO TO 230
            KM1 = K - 1
            KK = IK + K
            IKP1 = IK + K
            IF (KPVT(K) .LT. 0) GO TO 150
C
C              1 BY 1
C
               AP(KK) = (1.0E0,0.0E0)/AP(KK)
               IF (KM1 .LT. 1) GO TO 140
                  CALL CCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 130 J = 1, KM1
                     JK = IK + J
                     AP(JK) = CDOTU(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  130             CONTINUE
                  AP(KK) = AP(KK) + CDOTU(KM1,WORK,1,AP(IK+1),1)
  140          CONTINUE
               KSTEP = 1
            GO TO 190
  150       CONTINUE
C
C              2 BY 2
C
               KKP1 = IKP1 + K
               T = AP(KKP1)
               AK = AP(KK)/T
               AKP1 = AP(KKP1+1)/T
               AKKP1 = AP(KKP1)/T
               D = T*(AK*AKP1 - (1.0E0,0.0E0))
               AP(KK) = AKP1/D
               AP(KKP1+1) = AK/D
               AP(KKP1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 180
                  CALL CCOPY(KM1,AP(IKP1+1),1,WORK,1)
                  IJ = 0
                  DO 160 J = 1, KM1
                     JKP1 = IKP1 + J
                     AP(JKP1) = CDOTU(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                     IJ = IJ + J
  160             CONTINUE
                  AP(KKP1+1) = AP(KKP1+1)
     1                         + CDOTU(KM1,WORK,1,AP(IKP1+1),1)
                  AP(KKP1) = AP(KKP1)
     1                       + CDOTU(KM1,AP(IK+1),1,AP(IKP1+1),1)
                  CALL CCOPY(KM1,AP(IK+1),1,WORK,1)
                  IJ = 0
                  DO 170 J = 1, KM1
                     JK = IK + J
                     AP(JK) = CDOTU(J,AP(IJ+1),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                     IJ = IJ + J
  170             CONTINUE
                  AP(KK) = AP(KK) + CDOTU(KM1,WORK,1,AP(IK+1),1)
  180          CONTINUE
               KSTEP = 2
  190       CONTINUE
C
C           SWAP
C
            KS = ABS(KPVT(K))
            IF (KS .EQ. K) GO TO 220
               IKS = (KS*(KS - 1))/2
               CALL CSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
               KSJ = IK + KS
               DO 200 JB = KS, K
                  J = K + KS - JB
                  JK = IK + J
                  TEMP = AP(JK)
                  AP(JK) = AP(KSJ)
                  AP(KSJ) = TEMP
                  KSJ = KSJ - (J - 1)
  200          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 210
                  KSKP1 = IKP1 + KS
                  TEMP = AP(KSKP1)
                  AP(KSKP1) = AP(KKP1)
                  AP(KKP1) = TEMP
  210          CONTINUE
  220       CONTINUE
            IK = IK + K
            IF (KSTEP .EQ. 2) IK = IK + K + 1
            K = K + KSTEP
         GO TO 120
  230    CONTINUE
  240 CONTINUE
      RETURN
      END
