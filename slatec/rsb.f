*DECK RSB
      SUBROUTINE RSB (NM, N, MB, A, W, MATZ, Z, FV1, FV2, IERR)
C***BEGIN PROLOGUE  RSB
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a symmetric band matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A6
C***TYPE      SINGLE PRECISION (RSB-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC BAND matrix.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        MB is the half band width of the matrix, defined as the
C          number of adjacent diagonals, including the principal
C          diagonal, required to specify the non-zero portion of the
C          lower triangle of the matrix.  MB must be less than or
C          equal to N.  MB is an INTEGER variable.
C
C        A contains the lower triangle of the real symmetric band
C          matrix.  Its lowest subdiagonal is stored in the last
C          N+1-MB  positions of the first column, its next subdiagonal
C          in the last  N+2-MB  positions of the second column, further
C          subdiagonals similarly, and finally its principal diagonal
C          in the  N  positions of the last column.  Contents of storage
C          locations not part of the matrix are arbitrary.  A is a
C          two-dimensional REAL array, dimensioned A(NM,MB).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        A has been destroyed.
C
C        W contains the eigenvalues in ascending order.  W is a one-
C          dimensional REAL array, dimensioned W(N).
C
C        Z contains the eigenvectors if MATZ is not zero.  The
C          eigenvectors are orthonormal.  Z is a two-dimensional
C          REAL array, dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          12*N       if MB is either non-positive or greater than N,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues and eigenvectors, if requested,
C                     should be correct for indices 1, 2, ..., IERR-1.
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
C***ROUTINES CALLED  BANDR, TQL2, TQLRAT
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RSB
C
      INTEGER N,MB,NM,IERR,MATZ
      REAL A(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)
      LOGICAL TF
C
C***FIRST EXECUTABLE STATEMENT  RSB
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (MB .GT. 0) GO TO 10
      IERR = 12 * N
      GO TO 50
   10 IF (MB .LE. N) GO TO 15
      IERR = 12 * N
      GO TO 50
C
   15 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      TF = .FALSE.
      CALL  BANDR(NM,N,MB,A,W,FV1,FV2,TF,Z)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 TF = .TRUE.
      CALL  BANDR(NM,N,MB,A,W,FV1,FV1,TF,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
