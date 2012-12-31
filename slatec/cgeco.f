*DECK CGECO
      SUBROUTINE CGECO (A, LDA, N, IPVT, RCOND, Z)
C***BEGIN PROLOGUE  CGECO
C***PURPOSE  Factor a matrix using Gaussian elimination and estimate
C            the condition number of the matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1
C***TYPE      COMPLEX (SGECO-S, DGECO-D, CGECO-C)
C***KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGECO factors a complex matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, CGEFA is slightly faster.
C     To solve  A*X = B , follow CGECO By CGESL.
C     To Compute  INVERSE(A)*C , follow CGECO by CGESL.
C     To compute  DETERMINANT(A) , follow CGECO by CGEDI.
C     To compute  INVERSE(A) , follow CGECO by CGEDI.
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
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
C                in  A  and  B  of size  EPSILON  may cause
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
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CDOTC, CGEFA, CSSCAL, SCASUM
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGECO
      INTEGER LDA,N,IPVT(*)
      COMPLEX A(LDA,*),Z(*)
      REAL RCOND
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
      COMPLEX ZDUM,ZDUM1,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  CGECO
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = MAX(ANORM,SCASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL CGEFA(A,LDA,N,IPVT,INFO)
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
      DO 100 K = 1, N
         IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(A(K,K))) GO TO 30
            S = CABS1(A(K,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0E0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(A(K,K)) .EQ. 0.0E0) GO TO 40
            WK = WK/CONJG(A(K,K))
            WKM = WKM/CONJG(A(K,K))
         GO TO 50
   40    CONTINUE
            WK = (1.0E0,0.0E0)
            WKM = (1.0E0,0.0E0)
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))
               Z(J) = Z(J) + WK*CONJG(A(K,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*CONJG(A(K,J))
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
         IF (K .LT. N) Z(K) = Z(K) + CDOTC(N-K,A(K+1,K),1,Z(K+1),1)
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
         IF (K .LT. N) CALL CAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
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
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 150
            S = CABS1(A(K,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (CABS1(A(K,K)) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (CABS1(A(K,K)) .EQ. 0.0E0) Z(K) = (1.0E0,0.0E0)
         T = -Z(K)
         CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
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
