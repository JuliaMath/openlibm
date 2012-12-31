*DECK RSP
      SUBROUTINE RSP (NM, N, NV, A, W, MATZ, Z, FV1, FV2, IERR)
C***BEGIN PROLOGUE  RSP
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a real symmetric matrix packed into a one dimensional
C            array.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A1
C***TYPE      SINGLE PRECISION (RSP-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC PACKED matrix.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        NV is an INTEGER variable set equal to the dimension of the
C          array A as specified in the calling program.  NV must not
C          be less than  N*(N+1)/2.
C
C        A contains the lower triangle, stored row-wise, of the real
C          symmetric packed matrix.  A is a one-dimensional REAL
C          array, dimensioned A(NV).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        A has been destroyed.
C
C        W contains the eigenvalues in ascending order.  W is a
C          one-dimensional REAL array, dimensioned W(N).
C
C        Z contains the eigenvectors if MATZ is not zero.  The eigen-
C          vectors are orthonormal.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          20*N       if NV is less than N*(N+1)/2,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues and eigenvectors in the W and Z
C                     arrays should be correct for indices
C                     1, 2, ..., IERR-1.
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
C***ROUTINES CALLED  TQL2, TQLRAT, TRBAK3, TRED3
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RSP
C
      INTEGER I,J,N,NM,NV,IERR,MATZ
      REAL A(*),W(*),Z(NM,*),FV1(*),FV2(*)
C
C***FIRST EXECUTABLE STATEMENT  RSP
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10
      IERR = 20 * N
      GO TO 50
C
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)
      IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            Z(J,I) = 0.0E0
   30    CONTINUE
C
         Z(I,I) = 1.0E0
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  TRBAK3(NM,N,NV,A,N,Z)
   50 RETURN
      END
