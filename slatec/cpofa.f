*DECK CPOFA
      SUBROUTINE CPOFA (A, LDA, N, INFO)
C***BEGIN PROLOGUE  CPOFA
C***PURPOSE  Factor a complex Hermitian positive definite matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D1B
C***TYPE      COMPLEX (SPOFA-S, DPOFA-D, CPOFA-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
C             POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CPOFA factors a complex Hermitian positive definite matrix.
C
C     CPOFA is usually called by CPOCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for CPOCO) = (1 + 18/N)*(Time for CPOFA) .
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
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CDOTC
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CPOFA
      INTEGER LDA,N,INFO
      COMPLEX A(LDA,*)
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER J,JM1,K
C***FIRST EXECUTABLE STATEMENT  CPOFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - CDOTC(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + REAL(T*CONJG(T))
   10       CONTINUE
   20       CONTINUE
            S = REAL(A(J,J)) - S
            IF (S .LE. 0.0E0 .OR. AIMAG(A(J,J)) .NE. 0.0E0) GO TO 40
            A(J,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
