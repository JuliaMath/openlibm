*DECK SOSSOL
      SUBROUTINE SOSSOL (K, N, L, X, C, B, M)
C***BEGIN PROLOGUE  SOSSOL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SOS
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SOSSOL-S, DSOSSL-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     SOSSOL solves an upper triangular type of linear system by back
C     substitution.
C
C     The matrix C is upper trapezoidal and stored as a linear array by
C     rows. The equations have been normalized so that the diagonal
C     entries of C are understood to be unity. The off diagonal entries
C     and the elements of the constant right hand side vector B have
C     already been stored as the negatives of the corresponding equation
C     values.
C     with each call to SOSSOL a (K-1) by (K-1) triangular system is
C     resolved. For L greater than K, column L of C is included in the
C     right hand side vector.
C
C***SEE ALSO  SOS
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   801001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  SOSSOL
C
C
      DIMENSION X(*), C(*), B(*)
C
C***FIRST EXECUTABLE STATEMENT  SOSSOL
      NP1 = N + 1
      KM1 = K - 1
      LK = KM1
      IF (L .EQ. K) LK = K
      KN = M
C
C
      DO 40 KJ=1,KM1
        KMM1 = K - KJ
        KM = KMM1 + 1
        XMAX = 0.
        KN = KN - NP1 + KMM1
        IF (KM .GT. LK) GO TO 20
        JKM = KN
C
        DO 10 J=KM,LK
          JKM = JKM + 1
          XMAX = XMAX + C(JKM)*X(J)
   10   CONTINUE
C
   20   IF (L .LE. K) GO TO 30
        JKM = KN + L - KMM1
        XMAX = XMAX + C(JKM)*X(L)
   30   X(KMM1) = XMAX + B(KMM1)
   40 CONTINUE
C
      RETURN
      END
