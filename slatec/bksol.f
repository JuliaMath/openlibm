*DECK BKSOL
      SUBROUTINE BKSOL (N, A, X)
C***BEGIN PROLOGUE  BKSOL
C***SUBSIDIARY
C***PURPOSE  Subsidiary to BVSUP
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (BKSOL-S, DBKSOL-D)
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
C***SEE ALSO  BVSUP
C***ROUTINES CALLED  SDOT
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  BKSOL
C
      DIMENSION A(*),X(*)
C
C***FIRST EXECUTABLE STATEMENT  BKSOL
      M=(N*(N+1))/2
      X(N)=X(N)*A(M)
      IF (N .EQ. 1) GO TO 20
      NM1=N-1
      DO 10 K=1,NM1
      J=N-K
      M=M-K-1
   10 X(J)=X(J)*A(M) - SDOT(K,A(M+1),1,X(J+1),1)
C
   20 RETURN
      END
