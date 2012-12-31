*DECK CPOCO
      SUBROUTINE CPOCO (A, LDA, N, RCOND, Z, INFO)
C***BEGIN PROLOGUE  CPOCO
C***PURPOSE  Factor a complex Hermitian positive definite matrix
C            and estimate the condition number of the matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D1B
C***TYPE      COMPLEX (SPOCO-S, DPOCO-D, CPOCO-C)
C***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION, POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CPOCO factors a complex Hermitian positive definite matrix
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, CPOFA is slightly faster.
C     To solve  A*X = B , follow CPOCO by CPOSL.
C     To compute  INVERSE(A)*C , follow CPOCO by CPOSL.
C     To compute  DETERMINANT(A) , follow CPOCO by CPODI.
C     To compute  INVERSE(A) , follow CPOCO by CPODI.
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the Hermitian matrix to be factored.  Only the
C                diagonal and upper triangle are used.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix  R  so that  A =
C                CTRANS(R)*R where  CTRANS(R)  is the conjugate
C                transpose.  The strict lower triangle is unaltered.
C                If  INFO .NE. 0 , the factorization is not complete.
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
C                underflows.  If INFO .NE. 0 , RCOND is unchanged.
C
C        Z       COMPLEX(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                If  INFO .NE. 0 , Z  is unchanged.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CDOTC, CPOFA, CSSCAL, SCASUM
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CPOCO
      INTEGER LDA,N,INFO
      COMPLEX A(LDA,*),Z(*)
      REAL RCOND
C
      COMPLEX CDOTC,EK,T,WK,WKM
      REAL ANORM,S,SCASUM,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
      COMPLEX ZDUM,ZDUM2,CSIGN1
      REAL CABS1
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
      CSIGN1(ZDUM,ZDUM2) = CABS1(ZDUM)*(ZDUM2/CABS1(ZDUM2))
C
C     FIND NORM OF A USING ONLY UPPER HALF
C
C***FIRST EXECUTABLE STATEMENT  CPOCO
      DO 30 J = 1, N
         Z(J) = CMPLX(SCASUM(J,A(1,J),1),0.0E0)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = CMPLX(REAL(Z(I))+CABS1(A(I,J)),0.0E0)
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0E0
      DO 40 J = 1, N
         ANORM = MAX(ANORM,REAL(Z(J)))
   40 CONTINUE
C
C     FACTOR
C
      CALL CPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE CTRANS(R)*W = E
C
         EK = (1.0E0,0.0E0)
         DO 50 J = 1, N
            Z(J) = (0.0E0,0.0E0)
   50    CONTINUE
         DO 110 K = 1, N
            IF (CABS1(Z(K)) .NE. 0.0E0) EK = CSIGN1(EK,-Z(K))
            IF (CABS1(EK-Z(K)) .LE. REAL(A(K,K))) GO TO 60
               S = REAL(A(K,K))/CABS1(EK-Z(K))
               CALL CSSCAL(N,S,Z,1)
               EK = CMPLX(S,0.0E0)*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = CABS1(WK)
            SM = CABS1(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))
                  Z(J) = Z(J) + WK*CONJG(A(K,J))
                  S = S + CABS1(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*CONJG(A(K,J))
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
C        SOLVE R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(A(K,K))) GO TO 120
               S = REAL(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
C
         YNORM = 1.0E0
C
C        SOLVE CTRANS(R)*V = Y
C
         DO 150 K = 1, N
            Z(K) = Z(K) - CDOTC(K-1,A(1,K),1,Z(1),1)
            IF (CABS1(Z(K)) .LE. REAL(A(K,K))) GO TO 140
               S = REAL(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE R*Z = V
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (CABS1(Z(K)) .LE. REAL(A(K,K))) GO TO 160
               S = REAL(A(K,K))/CABS1(Z(K))
               CALL CSSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0E0/SCASUM(N,Z,1)
         CALL CSSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
  180 CONTINUE
      RETURN
      END
