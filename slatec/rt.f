*DECK RT
      SUBROUTINE RT (NM, N, A, W, MATZ, Z, FV1, IERR)
C***BEGIN PROLOGUE  RT
C***PURPOSE  Compute the eigenvalues and eigenvectors of a special real
C            tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5
C***TYPE      SINGLE PRECISION (RT-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of subroutines
C     from the eigensystem subroutine package (EISPACK) to find the
C     eigenvalues and eigenvectors (if desired) of a special REAL
C     TRIDIAGONAL matrix.  The property of the matrix required for use
C     of this subroutine is that the products of pairs of corresponding
C     off-diagonal elements be all non-negative.  If eigenvectors are
C     desired, no product can be zero unless both factors are zero.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        A contains the special real tridiagonal matrix in its first
C          three columns.  The subdiagonal elements are stored in the
C          last N-1 positions of the first column, the diagonal elements
C          in the second column, and the superdiagonal elements in the
C          first N-1 positions of the third column.  Elements A(1,1) and
C          A(N,3) are arbitrary.  A is a two-dimensional REAL array,
C          dimensioned A(NM,3).
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
C        Z contains the eigenvectors if MATZ is not zero.  The eigen-
C          vectors are not normalized.  Z is a two-dimensional REAL
C          array, dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          N+J        if A(J,1)*A(J-1,3) is negative,
C          2*N+J      if the product is zero with one factor non-zero,
C                     and MATZ is non-zero;
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues and eigenvectors in the W and Z
C                     arrays should be correct for indices
C                     1, 2, ..., IERR-1.
C
C        FV1 is a one-dimensional REAL array used for temporary storage,
C          dimensioned FV1(N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  FIGI, FIGI2, IMTQL1, IMTQL2
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RT
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,3),W(*),Z(NM,*),FV1(*)
C
C***FIRST EXECUTABLE STATEMENT  RT
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  FIGI(NM,N,A,W,FV1,FV1,IERR)
      IF (IERR .GT. 0) GO TO 50
      CALL  IMTQL1(N,W,FV1,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  FIGI2(NM,N,A,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  IMTQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
