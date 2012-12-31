*DECK CHPFA
      SUBROUTINE CHPFA (AP, N, KPVT, INFO)
C***BEGIN PROLOGUE  CHPFA
C***PURPOSE  Factor a complex Hermitian matrix stored in packed form by
C            elimination with symmetric pivoting.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D1A
C***TYPE      COMPLEX (SSPFA-S, DSPFA-D, CHPFA-C, DSPFA-C)
C***KEYWORDS  HERMITIAN, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
C             PACKED
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     CHPFA factors a complex Hermitian matrix stored in
C     packed form by elimination with symmetric pivoting.
C
C     To solve  A*X = B , follow CHPFA by CHPSL.
C     To compute  INVERSE(A)*C , follow CHPFA by CHPSL.
C     To compute  DETERMINANT(A) , follow CHPFA by CHPDI.
C     To compute  INERTIA(A) , follow CHPFA by CHPDI.
C     To compute  INVERSE(A) , follow CHPFA by CHPDI.
C
C     On Entry
C
C        AP      COMPLEX (N*(N+1)/2)
C                the packed form of a Hermitian matrix  A .  The
C                columns of the upper triangle are stored sequentially
C                in a one-dimensional array of length  N*(N+1)/2 .
C                See comments below for details.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     Output
C
C        AP      A block diagonal matrix and the multipliers which
C                were used to obtain it stored in packed form.
C                The factorization can be written  A = U*D*CTRANS(U)
C                where  U  is a product of permutation and unit
C                upper triangular matrices , CTRANS(U) is the
C                conjugate transpose of  U , and  D  is block diagonal
C                with 1 by 1 and 2 by 2 blocks.
C
C        KVPT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if the K-th pivot block is singular.  This is
C                     not an error condition for this subroutine,
C                     but it does indicate that CHPSL or CHPDI may
C                     divide by zero if called.
C
C     Packed Storage
C
C          The following program segment will pack the upper
C          triangle of a Hermitian matrix.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K)  = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSWAP, ICAMAX
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
C***END PROLOGUE  CHPFA
      INTEGER N,KPVT(*),INFO
      COMPLEX AP(*)
C
      COMPLEX AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER ICAMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CHPFA
C
C     INITIALIZE
C
C     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
C
      ALPHA = (1.0E0 + SQRT(17.0E0))/8.0E0
C
      INFO = 0
C
C     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
C
      K = N
      IK = (N*(N - 1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (CABS1(AP(1)) .EQ. 0.0E0) INFO = 1
            GO TO 200
   20    CONTINUE
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS
C        REQUIRED.
C
         KM1 = K - 1
         KK = IK + K
         ABSAKK = CABS1(AP(KK))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = ICAMAX(K-1,AP(IK+1),1)
         IMK = IK + IMAX
         COLMAX = CABS1(AP(IMK))
         IF (ABSAKK .LT. ALPHA*COLMAX) GO TO 30
            KSTEP = 1
            SWAP = .FALSE.
         GO TO 90
   30    CONTINUE
C
C           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C           ROW IMAX.
C
            ROWMAX = 0.0E0
            IMAXP1 = IMAX + 1
            IM = IMAX*(IMAX - 1)/2
            IMJ = IM + 2*IMAX
            DO 40 J = IMAXP1, K
               ROWMAX = MAX(ROWMAX,CABS1(AP(IMJ)))
               IMJ = IMJ + J
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ICAMAX(IMAX-1,AP(IM+1),1)
               JMIM = JMAX + IM
               ROWMAX = MAX(ROWMAX,CABS1(AP(JMIM)))
   50       CONTINUE
            IMIM = IMAX + IM
            IF (CABS1(AP(IMIM)) .LT. ALPHA*ROWMAX) GO TO 60
               KSTEP = 1
               SWAP = .TRUE.
            GO TO 80
   60       CONTINUE
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) GO TO 70
               KSTEP = 1
               SWAP = .FALSE.
            GO TO 80
   70       CONTINUE
               KSTEP = 2
               SWAP = IMAX .NE. KM1
   80       CONTINUE
   90    CONTINUE
         IF (MAX(ABSAKK,COLMAX) .NE. 0.0E0) GO TO 100
C
C           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
C
            KPVT(K) = K
            INFO = K
         GO TO 190
  100    CONTINUE
         IF (KSTEP .EQ. 2) GO TO 140
C
C           1 X 1 PIVOT BLOCK.
C
            IF (.NOT.SWAP) GO TO 120
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
               IMJ = IK + IMAX
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  JK = IK + J
                  T = CONJG(AP(JK))
                  AP(JK) = CONJG(AP(IMJ))
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            IJ = IK - (K - 1)
            DO 130 JJ = 1, KM1
               J = K - JJ
               JK = IK + J
               MULK = -AP(JK)/AP(KK)
               T = CONJG(MULK)
               CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
               IJJ = IJ + J
               AP(IJJ) = CMPLX(REAL(AP(IJJ)),0.0E0)
               AP(JK) = MULK
               IJ = IJ - (J - 1)
  130       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = K
            IF (SWAP) KPVT(K) = IMAX
         GO TO 190
  140    CONTINUE
C
C           2 X 2 PIVOT BLOCK.
C
            KM1K = IK + K - 1
            IKM1 = IK - (K - 1)
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
               IMJ = IKM1 + IMAX
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  JKM1 = IKM1 + J
                  T = CONJG(AP(JKM1))
                  AP(JKM1) = CONJG(AP(IMJ))
                  AP(IMJ) = T
                  IMJ = IMJ - (J - 1)
  150          CONTINUE
               T = AP(KM1K)
               AP(KM1K) = AP(IMK)
               AP(IMK) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = AP(KK)/AP(KM1K)
               KM1KM1 = IKM1 + K - 1
               AKM1 = AP(KM1KM1)/CONJG(AP(KM1K))
               DENOM = 1.0E0 - AK*AKM1
               IJ = IK - (K - 1) - (K - 2)
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  JK = IK + J
                  BK = AP(JK)/AP(KM1K)
                  JKM1 = IKM1 + J
                  BKM1 = AP(JKM1)/CONJG(AP(KM1K))
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = CONJG(MULK)
                  CALL CAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
                  T = CONJG(MULKM1)
                  CALL CAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
                  AP(JK) = MULK
                  AP(JKM1) = MULKM1
                  IJJ = IJ + J
                  AP(IJJ) = CMPLX(REAL(AP(IJJ)),0.0E0)
                  IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
