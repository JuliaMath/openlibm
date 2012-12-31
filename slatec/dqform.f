*DECK DQFORM
      SUBROUTINE DQFORM (M, N, Q, LDQ, WA)
C***BEGIN PROLOGUE  DQFORM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QFORM-S, DQFORM-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine proceeds from the computed QR factorization of
C     an M by N matrix A to accumulate the M by M orthogonal matrix
C     Q from its factored form.
C
C     The subroutine statement is
C
C       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A and the order of Q.
C
C       N is a positive integer input variable set to the number
C         of columns of A.
C
C       Q is an M by M array. On input the full lower trapezoid in
C         the first MIN(M,N) columns of Q contains the factored form.
C         On output Q has been accumulated into a square matrix.
C
C       LDQ is a positive integer input variable not less than M
C         which specifies the leading dimension of the array Q.
C
C       WA is a work array of length M.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQFORM
      INTEGER I, J, JM1, K, L, LDQ, M, MINMN, N, NP1
      DOUBLE PRECISION ONE, Q(LDQ,*), SUM, TEMP, WA(*), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
C
C***FIRST EXECUTABLE STATEMENT  DQFORM
      MINMN = MIN(M,N)
      IF (MINMN .LT. 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
C
      NP1 = N + 1
      IF (M .LT. NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
C
C     ACCUMULATE Q FROM ITS FACTORED FORM.
C
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (WA(K) .EQ. ZERO) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DQFORM.
C
      END
