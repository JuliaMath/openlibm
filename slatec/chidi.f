*DECK CHIDI
      SUBROUTINE CHIDI (A, LDA, N, KPVT, DET, INERT, WORK, JOB)
C***BEGIN PROLOGUE  CHIDI
C***PURPOSE  Compute the determinant, inertia and inverse of a complex
C            Hermitian matrix using the factors obtained from CHIFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D1A, D3D1A
C***TYPE      COMPLEX (SSIDI-S, DSISI-D, CHIDI-C, CSIDI-C)
C***KEYWORDS  DETERMINANT, HERMITIAN, INVERSE, LINEAR ALGEBRA, LINPACK,
C             MATRIX
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     CHIDI computes the determinant, inertia and inverse
C     of a complex Hermitian matrix using the factors from CHIFA.
C
C     On Entry
C
C        A       COMPLEX(LDA,N)
C                the output from CHIFA.
C
C        LDA     INTEGER
C                the leading dimension of the array A.
C
C        N       INTEGER
C                the order of the matrix A.
C
C        KVPT    INTEGER(N)
C                the pivot vector from CHIFA.
C
C        WORK    COMPLEX(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                JOB has the decimal expansion  ABC  where
C                   if  C .NE. 0, the inverse is computed,
C                   if  B .NE. 0, the determinant is computed,
C                   if  A .NE. 0, the inertia is computed.
C
C                For example, JOB = 111  gives all three.
C
C     On Return
C
C        Variables not requested by JOB are not used.
C
C        A      contains the upper triangle of the inverse of
C               the original matrix.  The strict lower triangle
C               is never referenced.
C
C        DET    REAL(2)
C               determinant of original matrix.
C               Determinant = DET(1) * 10.0**DET(2)
C               with 1.0 .LE. ABS(DET(1)) .LT. 10.0
C               or DET(1) = 0.0.
C
C        INERT  INTEGER(3)
C               the inertia of the original matrix.
C               INERT(1)  =  number of positive eigenvalues.
C               INERT(2)  =  number of negative eigenvalues.
C               INERT(3)  =  number of zero eigenvalues.
C
C     Error Condition
C
C        A division by zero may occur if the inverse is requested
C        and  CHICO  has set RCOND .EQ. 0.0
C        or  CHIFA  has set  INFO .NE. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CCOPY, CDOTC, CSWAP
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
C***END PROLOGUE  CHIDI
      INTEGER LDA,N,JOB
      COMPLEX A(LDA,*),WORK(*)
      REAL DET(2)
      INTEGER KPVT(*),INERT(3)
C
      COMPLEX AKKP1,CDOTC,TEMP
      REAL TEN,D,T,AK,AKP1
      INTEGER J,JB,K,KM1,KS,KSTEP
      LOGICAL NOINV,NODET,NOERT
C***FIRST EXECUTABLE STATEMENT  CHIDI
      NOINV = MOD(JOB,10) .EQ. 0
      NODET = MOD(JOB,100)/10 .EQ. 0
      NOERT = MOD(JOB,1000)/100 .EQ. 0
C
      IF (NODET .AND. NOERT) GO TO 140
         IF (NOERT) GO TO 10
            INERT(1) = 0
            INERT(2) = 0
            INERT(3) = 0
   10    CONTINUE
         IF (NODET) GO TO 20
            DET(1) = 1.0E0
            DET(2) = 0.0E0
            TEN = 10.0E0
   20    CONTINUE
         T = 0.0E0
         DO 130 K = 1, N
            D = REAL(A(K,K))
C
C           CHECK IF 1 BY 1
C
            IF (KPVT(K) .GT. 0) GO TO 50
C
C              2 BY 2 BLOCK
C              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S)
C                      (S  C)
C              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
C              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
C
               IF (T .NE. 0.0E0) GO TO 30
                  T = ABS(A(K,K+1))
                  D = (D/T)*REAL(A(K+1,K+1)) - T
               GO TO 40
   30          CONTINUE
                  D = T
                  T = 0.0E0
   40          CONTINUE
   50       CONTINUE
C
            IF (NOERT) GO TO 60
               IF (D .GT. 0.0E0) INERT(1) = INERT(1) + 1
               IF (D .LT. 0.0E0) INERT(2) = INERT(2) + 1
               IF (D .EQ. 0.0E0) INERT(3) = INERT(3) + 1
   60       CONTINUE
C
            IF (NODET) GO TO 120
               DET(1) = D*DET(1)
               IF (DET(1) .EQ. 0.0E0) GO TO 110
   70             IF (ABS(DET(1)) .GE. 1.0E0) GO TO 80
                     DET(1) = TEN*DET(1)
                     DET(2) = DET(2) - 1.0E0
                  GO TO 70
   80             CONTINUE
   90             IF (ABS(DET(1)) .LT. TEN) GO TO 100
                     DET(1) = DET(1)/TEN
                     DET(2) = DET(2) + 1.0E0
                  GO TO 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C     COMPUTE INVERSE(A)
C
      IF (NOINV) GO TO 270
         K = 1
  150    IF (K .GT. N) GO TO 260
            KM1 = K - 1
            IF (KPVT(K) .LT. 0) GO TO 180
C
C              1 BY 1
C
               A(K,K) = CMPLX(1.0E0/REAL(A(K,K)),0.0E0)
               IF (KM1 .LT. 1) GO TO 170
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 160 J = 1, KM1
                     A(J,K) = CDOTC(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
                  A(K,K) = A(K,K)
     1                     + CMPLX(REAL(CDOTC(KM1,WORK,1,A(1,K),1)),
     2                             0.0E0)
  170          CONTINUE
               KSTEP = 1
            GO TO 220
  180       CONTINUE
C
C              2 BY 2
C
               T = ABS(A(K,K+1))
               AK = REAL(A(K,K))/T
               AKP1 = REAL(A(K+1,K+1))/T
               AKKP1 = A(K,K+1)/T
               D = T*(AK*AKP1 - 1.0E0)
               A(K,K) = CMPLX(AKP1/D,0.0E0)
               A(K+1,K+1) = CMPLX(AK/D,0.0E0)
               A(K,K+1) = -AKKP1/D
               IF (KM1 .LT. 1) GO TO 210
                  CALL CCOPY(KM1,A(1,K+1),1,WORK,1)
                  DO 190 J = 1, KM1
                     A(J,K+1) = CDOTC(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
                  A(K+1,K+1) = A(K+1,K+1)
     1                         + CMPLX(REAL(CDOTC(KM1,WORK,1,A(1,K+1),
     2                                            1)),0.0E0)
                  A(K,K+1) = A(K,K+1) + CDOTC(KM1,A(1,K),1,A(1,K+1),1)
                  CALL CCOPY(KM1,A(1,K),1,WORK,1)
                  DO 200 J = 1, KM1
                     A(J,K) = CDOTC(J,A(1,J),1,WORK,1)
                     CALL CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
                  A(K,K) = A(K,K)
     1                     + CMPLX(REAL(CDOTC(KM1,WORK,1,A(1,K),1)),
     2                             0.0E0)
  210          CONTINUE
               KSTEP = 2
  220       CONTINUE
C
C           SWAP
C
            KS = ABS(KPVT(K))
            IF (KS .EQ. K) GO TO 250
               CALL CSWAP(KS,A(1,KS),1,A(1,K),1)
               DO 230 JB = KS, K
                  J = K + KS - JB
                  TEMP = CONJG(A(J,K))
                  A(J,K) = CONJG(A(KS,J))
                  A(KS,J) = TEMP
  230          CONTINUE
               IF (KSTEP .EQ. 1) GO TO 240
                  TEMP = A(KS,K+1)
                  A(KS,K+1) = A(K,K+1)
                  A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
            K = K + KSTEP
         GO TO 150
  260    CONTINUE
  270 CONTINUE
      RETURN
      END
