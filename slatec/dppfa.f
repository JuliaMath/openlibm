*DECK DPPFA
      SUBROUTINE DPPFA (AP, N, INFO)
C***BEGIN PROLOGUE  DPPFA
C***PURPOSE  Factor a real symmetric positive definite matrix stored in
C            packed form.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B1B
C***TYPE      DOUBLE PRECISION (SPPFA-S, DPPFA-D, CPPFA-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
C             POSITIVE DEFINITE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DPPFA factors a double precision symmetric positive definite
C     matrix stored in packed form.
C
C     DPPFA is usually called by DPPCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (time for DPPCO) = (1 + 18/N)*(time for DPPFA) .
C
C     On Entry
C
C        AP      DOUBLE PRECISION (N*(N+1)/2)
C                the packed form of a symmetric matrix  A .  The
C                columns of the upper triangle are stored sequentially
C                in a one-dimensional array of length  N*(N+1)/2 .
C                See comments below for details.
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        AP      an upper triangular matrix  R , stored in packed
C                form, so that  A = TRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  for normal return.
C                = K  if the leading minor of order  K  is not
C                     positive definite.
C
C
C     Packed Storage
C
C          The following program segment will pack the upper
C          triangle of a symmetric matrix.
C
C                K = 0
C                DO 20 J = 1, N
C                   DO 10 I = 1, J
C                      K = K + 1
C                      AP(K) = A(I,J)
C             10    CONTINUE
C             20 CONTINUE
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
C***END PROLOGUE  DPPFA
      INTEGER N,INFO
      DOUBLE PRECISION AP(*)
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER J,JJ,JM1,K,KJ,KK
C***FIRST EXECUTABLE STATEMENT  DPPFA
         JJ = 0
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            KJ = JJ
            KK = 0
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               KJ = KJ + 1
               T = AP(KJ) - DDOT(K-1,AP(KK+1),1,AP(JJ+1),1)
               KK = KK + K
               T = T/AP(KK)
               AP(KJ) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            JJ = JJ + J
            S = AP(JJ) - S
            IF (S .LE. 0.0D0) GO TO 40
            AP(JJ) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
