*DECK RST
      SUBROUTINE RST (NM, N, W, E, MATZ, Z, IERR)
C***BEGIN PROLOGUE  RST
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a real symmetric tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5
C***TYPE      SINGLE PRECISION (RST-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a REAL SYMMETRIC TRIDIAGONAL matrix.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        W contains the diagonal elements of the real symmetric
C          tridiagonal matrix.  W is a one-dimensional REAL array,
C          dimensioned W(N).
C
C        E contains the subdiagonal elements of the matrix in its last
C          N-1 positions.  E(1) is arbitrary.  E is a one-dimensional
C          REAL array, dimensioned E(N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        W contains the eigenvalues in ascending order.
C
C        Z contains the eigenvectors if MATZ is not zero.  The eigen-
C          vectors are orthonormal.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues and eigenvectors in the W and Z
C                     arrays should be correct for indices
C                     1, 2, ..., IERR-1.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  IMTQL1, IMTQL2
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RST
C
      INTEGER I,J,N,NM,IERR,MATZ
      REAL W(*),E(*),Z(NM,*)
C
C***FIRST EXECUTABLE STATEMENT  RST
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
      CALL  IMTQL1(N,W,E,IERR)
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
      CALL  IMTQL2(NM,N,W,E,Z,IERR)
   50 RETURN
      END
