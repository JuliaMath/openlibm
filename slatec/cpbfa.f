*DECK CPBFA
      SUBROUTINE CPBFA (ABD, LDA, N, M, INFO)
C***BEGIN PROLOGUE  CPBFA
C***PURPOSE  Factor a complex Hermitian positive definite matrix stored
C            in band form.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2D2
C***TYPE      COMPLEX (SPBFA-S, DPBFA-D, CPBFA-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
C             POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CPBFA factors a complex Hermitian positive definite matrix
C     stored in band form.
C
C     CPBFA is usually called by CPBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     COMPLEX(LDA, N)
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
C                form, so that  A = CTRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  if the leading minor of order  K  is not
C                     positive definite.
C
C     Band Storage
C
C           If  A  is a Hermitian positive definite band matrix,
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
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CDOTC
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CPBFA
      INTEGER LDA,N,M,INFO
      COMPLEX ABD(LDA,*)
C
      COMPLEX CDOTC,T
      REAL S
      INTEGER IK,J,JK,K,MU
C***FIRST EXECUTABLE STATEMENT  CPBFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0E0
            IK = M + 1
            JK = MAX(J-M,1)
            MU = MAX(M+2-J,1)
            IF (M .LT. MU) GO TO 20
            DO 10 K = MU, M
               T = ABD(K,J) - CDOTC(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
               T = T/ABD(M+1,JK)
               ABD(K,J) = T
               S = S + REAL(T*CONJG(T))
               IK = IK - 1
               JK = JK + 1
   10       CONTINUE
   20       CONTINUE
            S = REAL(ABD(M+1,J)) - S
            IF (S .LE. 0.0E0 .OR. AIMAG(ABD(M+1,J)) .NE. 0.0E0)
     1         GO TO 40
            ABD(M+1,J) = CMPLX(SQRT(S),0.0E0)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
