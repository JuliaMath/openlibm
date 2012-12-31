*DECK RSG
      SUBROUTINE RSG (NM, N, A, B, W, MATZ, Z, FV1, FV2, IERR)
C***BEGIN PROLOGUE  RSG
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a symmetric generalized eigenproblem.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4B1
C***TYPE      SINGLE PRECISION (RSG-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     To find the eigenvalues and eigenvectors (if desired)
C     for the REAL SYMMETRIC generalized eigenproblem  Ax = (LAMBDA)Bx.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A, B, and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrices A and B.  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        A contains a real symmetric matrix.  A is a two-dimensional
C          REAL array, dimensioned A(NM,N).
C
C        B contains a positive definite real symmetric matrix.  B is a
C          two-dimensional REAL array, dimensioned B(NM,N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        W contains the eigenvalues in ascending order.  W is a
C          one-dimensional REAL array, dimensioned W(N).
C
C        Z contains the eigenvectors if MATZ is not zero.  Z is a
C          two-dimensional REAL array, dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          7*N+1      if B is not positive definite,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues should be correct for indices
C                     1, 2, ..., IERR-1, but no eigenvectors are
C                     computed.
C
C        FV1 and FV2 are one-dimensional REAL arrays used for temporary
C          storage, dimensioned FV1(N) and FV2(N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  REBAK, REDUC, TQL2, TQLRAT, TRED1, TRED2
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RSG
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,*),B(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)
C
C***FIRST EXECUTABLE STATEMENT  RSG
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 CALL  REDUC(NM,N,A,B,FV2,IERR)
      IF (IERR .NE. 0) GO TO 50
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  REBAK(NM,N,B,FV2,N,Z)
   50 RETURN
      END
