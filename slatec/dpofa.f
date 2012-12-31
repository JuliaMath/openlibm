*DECK DPOFA
      SUBROUTINE DPOFA (A, LDA, N, INFO)
C***BEGIN PROLOGUE  DPOFA
C***PURPOSE  Factor a real symmetric positive definite matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1B
C***TYPE      DOUBLE PRECISION (SPOFA-S, DPOFA-D, CPOFA-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
C             POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DPOFA factors a double precision symmetric positive definite
C     matrix.
C
C     DPOFA is usually called by DPOCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (time for DPOCO) = (1 + 18/N)*(time for DPOFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the symmetric matrix to be factored.  Only the
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
C        A       an upper triangular matrix  R  so that  A = TRANS(R)*R
C                where  TRANS(R)  is the transpose.
C                The strict lower triangle is unaltered.
C                If  INFO .NE. 0 , the factorization is not complete.
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  signals an error condition.  The leading minor
C                     of order  K  is not positive definite.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPOFA
      INTEGER LDA,N,INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER J,JM1,K
C***FIRST EXECUTABLE STATEMENT  DPOFA
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
            IF (S .LE. 0.0D0) GO TO 40
            A(J,J) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
