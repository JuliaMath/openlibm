*DECK CPODI
      SUBROUTINE CPODI (A, LDA, N, DET, JOB)
C***BEGIN PROLOGUE  CPODI
C***PURPOSE  Compute the determinant and inverse of a certain complex
C            Hermitian positive definite matrix using the factors
C            computed by CPOCO, CPOFA, or CQRDC.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D1B, D3D1B
C***TYPE      COMPLEX (SPODI-S, DPODI-D, CPODI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
C             POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CPODI computes the determinant and inverse of a certain
C     complex Hermitian positive definite matrix (see below)
C     using the factors computed by CPOCO, CPOFA or CQRDC.
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the output  A  from CPOCO or CPOFA
C                or the output  X  from CQRDC.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       If CPOCO or CPOFA was used to factor  A  then
C                CPODI produces the upper half of INVERSE(A) .
C                If CQRDC was used to decompose  X  then
C                CPODI produces the upper half of INVERSE(CTRANS(X)*X)
C                where CTRANS(X) is the conjugate transpose.
C                Elements of  A  below the diagonal are unchanged.
C                If the units digit of JOB is zero,  A  is unchanged.
C
C        DET     REAL(2)
C                determinant of  A  or of  CTRANS(X)*X  if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DET(1) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        a division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if CPOCO or CPOFA has set INFO .EQ. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSCAL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CPODI
      INTEGER LDA,N,JOB
      COMPLEX A(LDA,*)
      REAL DET(2)
C
      COMPLEX T
      REAL S
      INTEGER I,J,JM1,K,KP1
C***FIRST EXECUTABLE STATEMENT  CPODI
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0E0
         DET(2) = 0.0E0
         S = 10.0E0
         DO 50 I = 1, N
            DET(1) = REAL(A(I,I))**2*DET(1)
            IF (DET(1) .EQ. 0.0E0) GO TO 60
   10       IF (DET(1) .GE. 1.0E0) GO TO 20
               DET(1) = S*DET(1)
               DET(2) = DET(2) - 1.0E0
            GO TO 10
   20       CONTINUE
   30       IF (DET(1) .LT. S) GO TO 40
               DET(1) = DET(1)/S
               DET(2) = DET(2) + 1.0E0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(R)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 140
         DO 100 K = 1, N
            A(K,K) = (1.0E0,0.0E0)/A(K,K)
            T = -A(K,K)
            CALL CSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = (0.0E0,0.0E0)
               CALL CAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM  INVERSE(R) * CTRANS(INVERSE(R))
C
         DO 130 J = 1, N
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 120
            DO 110 K = 1, JM1
               T = CONJG(A(K,J))
               CALL CAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
            T = CONJG(A(J,J))
            CALL CSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
      RETURN
      END
