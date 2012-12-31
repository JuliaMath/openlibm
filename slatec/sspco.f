*DECK SSPCO
      SUBROUTINE SSPCO (AP, N, KPVT, RCOND, Z)
C***BEGIN PROLOGUE  SSPCO
C***PURPOSE  Factor a real symmetric matrix stored in packed form
C            by elimination with symmetric pivoting and estimate the
C            condition number of the matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1A
C***TYPE      SINGLE PRECISION (SSPCO-S, DSPCO-D, CHPCO-C, CSPCO-C)
C***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION, PACKED, SYMMETRIC
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     SSPCO factors a real symmetric matrix stored in packed
C     form by elimination with symmetric pivoting and estimates
C     the condition of the matrix.
C
C     If  RCOND  is not needed, SSPFA is slightly faster.
C     To solve  A*X = B , follow SSPCO by SSPSL.
C     To compute  INVERSE(A)*C , follow SSPCO by SSPSL.
C     To compute  INVERSE(A) , follow SSPCO by SSPDI.
C     To compute  DETERMINANT(A) , follow SSPCO by SSPDI.
C     To compute  INERTIA(A), follow SSPCO by SSPDI.
C
C     On Entry
C
C        AP      REAL (N*(N+1)/2)
C                the packed form of a symmetric matrix  A .  The
C                columns of the upper triangle are stored sequentially
C                in a one-dimensional array of length  N*(N+1)/2 .
C                See comments below for details.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     Output
C
C        AP      a block diagonal matrix and the multipliers which
C                were used to obtain it stored in packed form.
C                The factorization can be written  A = U*D*TRANS(U)
C                where  U  is a product of permutation and unit
C                upper triangular matrices , TRANS(U) is the
C                transpose of  U , and  D  is block diagonal
C                with 1 by 1 and 2 by 2 blocks.
C
C        KPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       REAL(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     Packed Storage
C
C          The following program segment will pack the upper
C          triangle of a symmetric matrix.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SASUM, SAXPY, SDOT, SSCAL, SSPFA
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
C***END PROLOGUE  SSPCO
      INTEGER N,KPVT(*)
      REAL AP(*),Z(*)
      REAL RCOND
C
      REAL AK,AKM1,BK,BKM1,SDOT,DENOM,EK,T
      REAL ANORM,S,SASUM,YNORM
      INTEGER I,IJ,IK,IKM1,IKP1,INFO,J,JM1,J1
      INTEGER K,KK,KM1K,KM1KM1,KP,KPS,KS
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
C***FIRST EXECUTABLE STATEMENT  SSPCO
      J1 = 1
      DO 30 J = 1, N
         Z(J) = SASUM(J,AP(J1),1)
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + ABS(AP(IJ))
            IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
C
C     FACTOR
C
      CALL SSPFA(AP,N,KPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE U*D*W = E
C
      EK = 1.0E0
      DO 50 J = 1, N
         Z(J) = 0.0E0
   50 CONTINUE
      K = N
      IK = (N*(N - 1))/2
   60 IF (K .EQ. 0) GO TO 120
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         KP = ABS(KPVT(K))
         KPS = K + 1 - KS
         IF (KP .EQ. KPS) GO TO 70
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
   70    CONTINUE
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL SAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
         IF (KS .EQ. 1) GO TO 80
            IF (Z(K-1) .NE. 0.0E0) EK = SIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL SAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
   80    CONTINUE
         IF (KS .EQ. 2) GO TO 100
            IF (ABS(Z(K)) .LE. ABS(AP(KK))) GO TO 90
               S = ABS(AP(KK))/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
               EK = S*EK
   90       CONTINUE
            IF (AP(KK) .NE. 0.0E0) Z(K) = Z(K)/AP(KK)
            IF (AP(KK) .EQ. 0.0E0) Z(K) = 1.0E0
         GO TO 110
  100    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  110    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 60
  120 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
C     SOLVE TRANS(U)*Y = W
C
      K = 1
      IK = 0
  130 IF (K .GT. N) GO TO 160
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 150
            Z(K) = Z(K) + SDOT(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     1         Z(K+1) = Z(K+1) + SDOT(K-1,AP(IKP1+1),1,Z(1),1)
            KP = ABS(KPVT(K))
            IF (KP .EQ. K) GO TO 140
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  140       CONTINUE
  150    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 130
  160 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE U*D*V = Y
C
      K = N
      IK = N*(N - 1)/2
  170 IF (K .EQ. 0) GO TO 230
         KK = IK + K
         IKM1 = IK - (K - 1)
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. KS) GO TO 190
            KP = ABS(KPVT(K))
            KPS = K + 1 - KS
            IF (KP .EQ. KPS) GO TO 180
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
  180       CONTINUE
            CALL SAXPY(K-KS,Z(K),AP(IK+1),1,Z(1),1)
            IF (KS .EQ. 2) CALL SAXPY(K-KS,Z(K-1),AP(IKM1+1),1,Z(1),1)
  190    CONTINUE
         IF (KS .EQ. 2) GO TO 210
            IF (ABS(Z(K)) .LE. ABS(AP(KK))) GO TO 200
               S = ABS(AP(KK))/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
               YNORM = S*YNORM
  200       CONTINUE
            IF (AP(KK) .NE. 0.0E0) Z(K) = Z(K)/AP(KK)
            IF (AP(KK) .EQ. 0.0E0) Z(K) = 1.0E0
         GO TO 220
  210    CONTINUE
            KM1K = IK + K - 1
            KM1KM1 = IKM1 + K - 1
            AK = AP(KK)/AP(KM1K)
            AKM1 = AP(KM1KM1)/AP(KM1K)
            BK = Z(K)/AP(KM1K)
            BKM1 = Z(K-1)/AP(KM1K)
            DENOM = AK*AKM1 - 1.0E0
            Z(K) = (AKM1*BK - BKM1)/DENOM
            Z(K-1) = (AK*BKM1 - BK)/DENOM
  220    CONTINUE
         K = K - KS
         IK = IK - K
         IF (KS .EQ. 2) IK = IK - (K + 1)
      GO TO 170
  230 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE TRANS(U)*Z = V
C
      K = 1
      IK = 0
  240 IF (K .GT. N) GO TO 270
         KS = 1
         IF (KPVT(K) .LT. 0) KS = 2
         IF (K .EQ. 1) GO TO 260
            Z(K) = Z(K) + SDOT(K-1,AP(IK+1),1,Z(1),1)
            IKP1 = IK + K
            IF (KS .EQ. 2)
     1         Z(K+1) = Z(K+1) + SDOT(K-1,AP(IKP1+1),1,Z(1),1)
            KP = ABS(KPVT(K))
            IF (KP .EQ. K) GO TO 250
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
  250       CONTINUE
  260    CONTINUE
         IK = IK + K
         IF (KS .EQ. 2) IK = IK + (K + 1)
         K = K + KS
      GO TO 240
  270 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
