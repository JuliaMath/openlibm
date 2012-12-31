*DECK CHIFA
      SUBROUTINE CHIFA (A, LDA, N, KPVT, INFO)
C***BEGIN PROLOGUE  CHIFA
C***PURPOSE  Factor a complex Hermitian matrix by elimination
C            (symmetric pivoting).
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D1A
C***TYPE      COMPLEX (SSIFA-S, DSIFA-D, CHIFA-C, CSIFA-C)
C***KEYWORDS  HERMITIAN, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
C***AUTHOR  Bunch, J., (UCSD)
C***DESCRIPTION
C
C     CHIFA factors a complex Hermitian matrix by elimination
C     with symmetric pivoting.
C
C     To solve  A*X = B , follow CHIFA by CHISL.
C     To compute  INVERSE(A)*C , follow CHIFA by CHISL.
C     To compute  DETERMINANT(A) , follow CHIFA by CHIDI.
C     To compute  INERTIA(A) , follow CHIFA by CHIDI.
C     To compute  INVERSE(A) , follow CHIFA by CHIDI.
C
C     On Entry
C
C        A       COMPLEX(LDA,N)
C                the Hermitian matrix to be factored.
C                Only the diagonal and upper triangle are used.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       a block diagonal matrix and the multipliers which
C                were used to obtain it.
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
C                     but it does indicate that CHISL or CHIDI may
C                     divide by zero if called.
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
C***END PROLOGUE  CHIFA
      INTEGER LDA,N,KPVT(*),INFO
      COMPLEX A(LDA,*)
C
      COMPLEX AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
      REAL ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,ICAMAX
      LOGICAL SWAP
      COMPLEX ZDUM
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CHIFA
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
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
         IF (K .EQ. 0) GO TO 200
         IF (K .GT. 1) GO TO 20
            KPVT(1) = 1
            IF (CABS1(A(1,1)) .EQ. 0.0E0) INFO = 1
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
         ABSAKK = CABS1(A(K,K))
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
C        COLUMN K.
C
         IMAX = ICAMAX(K-1,A(1,K),1)
         COLMAX = CABS1(A(IMAX,K))
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
            DO 40 J = IMAXP1, K
               ROWMAX = MAX(ROWMAX,CABS1(A(IMAX,J)))
   40       CONTINUE
            IF (IMAX .EQ. 1) GO TO 50
               JMAX = ICAMAX(IMAX-1,A(1,IMAX),1)
               ROWMAX = MAX(ROWMAX,CABS1(A(JMAX,IMAX)))
   50       CONTINUE
            IF (CABS1(A(IMAX,IMAX)) .LT. ALPHA*ROWMAX) GO TO 60
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
               CALL CSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
               DO 110 JJ = IMAX, K
                  J = K + IMAX - JJ
                  T = CONJG(A(J,K))
                  A(J,K) = CONJG(A(IMAX,J))
                  A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            DO 130 JJ = 1, KM1
               J = K - JJ
               MULK = -A(J,K)/A(K,K)
               T = CONJG(MULK)
               CALL CAXPY(J,T,A(1,K),1,A(1,J),1)
               A(J,J) = CMPLX(REAL(A(J,J)),0.0E0)
               A(J,K) = MULK
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
            IF (.NOT.SWAP) GO TO 160
C
C              PERFORM AN INTERCHANGE.
C
               CALL CSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
               DO 150 JJ = IMAX, KM1
                  J = KM1 + IMAX - JJ
                  T = CONJG(A(J,K-1))
                  A(J,K-1) = CONJG(A(IMAX,J))
                  A(IMAX,J) = T
  150          CONTINUE
               T = A(K-1,K)
               A(K-1,K) = A(IMAX,K)
               A(IMAX,K) = T
  160       CONTINUE
C
C           PERFORM THE ELIMINATION.
C
            KM2 = K - 2
            IF (KM2 .EQ. 0) GO TO 180
               AK = A(K,K)/A(K-1,K)
               AKM1 = A(K-1,K-1)/CONJG(A(K-1,K))
               DENOM = 1.0E0 - AK*AKM1
               DO 170 JJ = 1, KM2
                  J = KM1 - JJ
                  BK = A(J,K)/A(K-1,K)
                  BKM1 = A(J,K-1)/CONJG(A(K-1,K))
                  MULK = (AKM1*BK - BKM1)/DENOM
                  MULKM1 = (AK*BKM1 - BK)/DENOM
                  T = CONJG(MULK)
                  CALL CAXPY(J,T,A(1,K),1,A(1,J),1)
                  T = CONJG(MULKM1)
                  CALL CAXPY(J,T,A(1,K-1),1,A(1,J),1)
                  A(J,K) = MULK
                  A(J,K-1) = MULKM1
                  A(J,J) = CMPLX(REAL(A(J,J)),0.0E0)
  170          CONTINUE
  180       CONTINUE
C
C           SET THE PIVOT ARRAY.
C
            KPVT(K) = 1 - K
            IF (SWAP) KPVT(K) = -IMAX
            KPVT(K-1) = KPVT(K)
  190    CONTINUE
         K = K - KSTEP
      GO TO 10
  200 CONTINUE
      RETURN
      END
