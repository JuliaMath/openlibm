*DECK SPBCO
      SUBROUTINE SPBCO (ABD, LDA, N, M, RCOND, Z, INFO)
C***BEGIN PROLOGUE  SPBCO
C***PURPOSE  Factor a real symmetric positive definite matrix stored in
C            band form and estimate the condition number of the matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B2
C***TYPE      SINGLE PRECISION (SPBCO-S, DPBCO-D, CPBCO-C)
C***KEYWORDS  BANDED, CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION, POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     SPBCO factors a real symmetric positive definite matrix
C     stored in band form and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, SPBFA is slightly faster.
C     To solve  A*X = B , follow SPBCO by SPBSL.
C     To compute  INVERSE(A)*C , follow SPBCO by SPBSL.
C     To compute  DETERMINANT(A) , follow SPBCO by SPBDI.
C
C     On Entry
C
C        ABD     REAL(LDA, N)
C                the matrix to be factored.  The columns of the upper
C                triangle are stored in the columns of ABD and the
C                diagonals of the upper triangle are stored in the
C                rows of ABD .  See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. M + 1 .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        M       INTEGER
C                the number of diagonals above the main diagonal.
C                0 .LE. M .LT. N .
C
C     On Return
C
C        ABD     an upper triangular matrix  R , stored in band
C                form, so that  A = TRANS(R)*R .
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
C        Z       REAL(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is singular to working precision, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C                If  INFO .NE. 0 , Z  is unchanged.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C     Band Storage
C
C           If  A  is a symmetric positive definite band matrix,
C           the following program segment will set up the input.
C
C                   M = (band width above diagonal)
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses  M + 1  rows of  A , except for the  M by M
C           upper left triangle, which is ignored.
C
C     Example:  If the original matrix is
C
C           11 12 13  0  0  0
C           12 22 23 24  0  0
C           13 23 33 34 35  0
C            0 24 34 44 45 46
C            0  0 35 45 55 56
C            0  0  0 46 56 66
C
C     then  N = 6 , M = 2  and  ABD  should contain
C
C            *  * 13 24 35 46
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  SASUM, SAXPY, SDOT, SPBFA, SSCAL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SPBCO
      INTEGER LDA,N,M,INFO
      REAL ABD(LDA,*),Z(*)
      REAL RCOND
C
      REAL SDOT,EK,T,WK,WKM
      REAL ANORM,S,SASUM,SM,YNORM
      INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU
C
C     FIND NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  SPBCO
      DO 30 J = 1, N
         L = MIN(J,M+1)
         MU = MAX(M+2-J,1)
         Z(J) = SASUM(L,ABD(MU,J),1)
         K = J - L
         IF (M .LT. MU) GO TO 20
         DO 10 I = MU, M
            K = K + 1
            Z(K) = Z(K) + ABS(ABD(I,J))
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
      CALL SPBFA(ABD,LDA,N,M,INFO)
      IF (INFO .NE. 0) GO TO 180
C
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C        SOLVE TRANS(R)*W = E
C
         EK = 1.0E0
         DO 50 J = 1, N
            Z(J) = 0.0E0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
            IF (ABS(EK-Z(K)) .LE. ABD(M+1,K)) GO TO 60
               S = ABD(M+1,K)/ABS(EK-Z(K))
               CALL SSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = ABS(WK)
            SM = ABS(WKM)
            WK = WK/ABD(M+1,K)
            WKM = WKM/ABD(M+1,K)
            KP1 = K + 1
            J2 = MIN(K+M,N)
            I = M + 1
            IF (KP1 .GT. J2) GO TO 100
               DO 70 J = KP1, J2
                  I = I - 1
                  SM = SM + ABS(Z(J)+WKM*ABD(I,J))
                  Z(J) = Z(J) + WK*ABD(I,J)
                  S = S + ABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  I = M + 1
                  DO 80 J = KP1, J2
                     I = I - 1
                     Z(J) = Z(J) + T*ABD(I,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
C
C        SOLVE  R*Y = W
C
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. ABD(M+1,K)) GO TO 120
               S = ABD(M+1,K)/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
  120       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL SAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  130    CONTINUE
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
C
         YNORM = 1.0E0
C
C        SOLVE TRANS(R)*V = Y
C
         DO 150 K = 1, N
            LM = MIN(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            Z(K) = Z(K) - SDOT(LM,ABD(LA,K),1,Z(LB),1)
            IF (ABS(Z(K)) .LE. ABD(M+1,K)) GO TO 140
               S = ABD(M+1,K)/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
  150    CONTINUE
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
C        SOLVE  R*Z = W
C
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. ABD(M+1,K)) GO TO 160
               S = ABD(M+1,K)/ABS(Z(K))
               CALL SSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/ABD(M+1,K)
            LM = MIN(K-1,M)
            LA = M + 1 - LM
            LB = K - LM
            T = -Z(K)
            CALL SAXPY(LM,T,ABD(LA,K),1,Z(LB),1)
  170    CONTINUE
C        MAKE ZNORM = 1.0
         S = 1.0E0/SASUM(N,Z,1)
         CALL SSCAL(N,S,Z,1)
         YNORM = S*YNORM
C
         IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
  180 CONTINUE
      RETURN
      END
