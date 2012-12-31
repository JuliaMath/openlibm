*DECK CGBCO
      SUBROUTINE CGBCO (ABD, LDA, N, ML, MU, IPVT, RCOND, Z)
C***BEGIN PROLOGUE  CGBCO
C***PURPOSE  Factor a band matrix by Gaussian elimination and
C            estimate the condition number of the matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C2
C***TYPE      COMPLEX (SGBCO-S, DGBCO-D, CGBCO-C)
C***KEYWORDS  BANDED, CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGBCO factors a complex band matrix by Gaussian
C     elimination and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, CGBFA is slightly faster.
C     To solve  A*X = B , follow CGBCO by CGBSL.
C     To compute  INVERSE(A)*C , follow CGBCO by CGBSL.
C     To compute  DETERMINANT(A) , follow CGBCO by CGBDI.
C
C     On Entry
C
C        ABD     COMPLEX(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT. N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT. N .
C                More efficient if  ML .LE. MU .
C
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  And  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       COMPLEX(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     Band Storage
C
C           if  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-MU)
C                      I2 = MIN(N, J+Ml)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C     Example:  If the original matrix is
C
C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66
C
C      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABD should contain
C
C            *  *  *  +  +  +  , * = not used
C            *  * 13 24 35 46  , + = used for pivoting
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CDOTC, CGBFA, CSSCAL, SCASUM
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGBCO
      INTEGER LDA,N,ML,MU,IPVT(*)
      COMPLEX ABD(LDA,*),Z(*)
      REAL RCOND
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  CGBCO
      ANORM = 0.0E0
      L = ML + 1
      IS = L + MU
      DO 10 J = 1, N
         ANORM = MAX(ANORM,SCASUM(L,ABD(IS,J),1))
         IF (IS .GT. ML + 1) IS = IS - 1
         IF (J .LE. MU) L = L + 1
         IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C
C     FACTOR
C
      CALL CGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(U)*W = E
C
      EK = (1.0E0,0.0E0)
      DO 20 J = 1, N
         Z(J) = (0.0E0,0.0E0)
   20 CONTINUE
      M = ML + MU + 1
      JU = 0
      DO 100 K = 1, N
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(ABD(M,K))) GO TO 30
            S = CABS1(ABD(M,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(ABD(M,K)) .EQ. 0.0E0) GO TO 40
            WK = WK/CONJG(ABD(M,K))
            WKM = WKM/CONJG(ABD(M,K))
         GO TO 50
   40    CONTINUE
            WK = (1.0E0,0.0E0)
            WKM = (1.0E0,0.0E0)
   50    CONTINUE
         KP1 = K + 1
         JU = MIN(MAX(JU,MU+IPVT(K)),N)
         MM = M
         IF (KP1 .GT. JU) GO TO 90
            DO 60 J = KP1, JU
               MM = MM - 1
               SM = SM + CABS1(Z(J)+WKM*CONJG(ABD(MM,J)))
               Z(J) = Z(J) + WK*CONJG(ABD(MM,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               MM = M
               DO 70 J = KP1, JU
                  MM = MM - 1
                  Z(J) = Z(J) + T*CONJG(ABD(MM,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         LM = MIN(ML,N-K)
         IF (K .LT. N) Z(K) = Z(K) + CDOTC(LM,ABD(M+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM = MIN(ML,N-K)
         IF (K .LT. N) CALL CAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = W
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (CABS1(Z(K)) .LE. CABS1(ABD(M,K))) GO TO 150
            S = CABS1(ABD(M,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (CABS1(ABD(M,K)) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
         IF (CABS1(ABD(M,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         LM = MIN(K,M) - 1
         LA = M - LM
         LZ = K - LM
         T = -Z(K)
         CALL CAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
