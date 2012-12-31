*DECK DBKSOL
      SUBROUTINE DBKSOL (N, A, X)
C***BEGIN PROLOGUE  DBKSOL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BKSOL-S, DBKSOL-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C     Solution of an upper triangular linear system by
C     back-substitution
C
C     The matrix A is assumed to be stored in a linear
C     array proceeding in a row-wise manner. The
C     vector X contains the given constant vector on input
C     and contains the solution on return.
C     The actual diagonal of A is unity while a diagonal
C     scaling matrix is stored there.
C **********************************************************************
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DBKSOL
C
      DOUBLE PRECISION DDOT
      INTEGER J, K, M, N, NM1
      DOUBLE PRECISION A(*), X(*)
C
C***FIRST EXECUTABLE STATEMENT  DBKSOL
      M = (N*(N + 1))/2
      X(N) = X(N)*A(M)
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 20
      DO 10 K = 1, NM1
         J = N - K
         M = M - K - 1
         X(J) = X(J)*A(M) - DDOT(K,A(M+1),1,X(J+1),1)
   10 CONTINUE
   20 CONTINUE
C
      RETURN
      END
