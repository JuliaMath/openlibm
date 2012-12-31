*DECK BDIFF
      SUBROUTINE BDIFF (L, V)
C***BEGIN PROLOGUE  BDIFF
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BSKIN
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BDIFF-S, DBDIFF-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     BDIFF computes the sum of B(L,K)*V(K)*(-1)**K where B(L,K)
C     are the binomial coefficients.  Truncated sums are computed by
C     setting last part of the V vector to zero. On return, the binomial
C     sum is in V(L).
C
C***SEE ALSO  BSKIN
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  BDIFF
      INTEGER I, J, K, L
      REAL V
      DIMENSION V(*)
C***FIRST EXECUTABLE STATEMENT  BDIFF
      IF (L.EQ.1) RETURN
      DO 20 J=2,L
        K = L
        DO 10 I=J,L
          V(K) = V(K-1) - V(K)
          K = K - 1
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
