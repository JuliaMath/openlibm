*DECK BAKVEC
      SUBROUTINE BAKVEC (NM, N, T, E, M, Z, IERR)
C***BEGIN PROLOGUE  BAKVEC
C***PURPOSE  Form the eigenvectors of a certain real non-symmetric
C            tridiagonal matrix from a symmetric tridiagonal matrix
C            output from FIGI.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (BAKVEC-S)
C***KEYWORDS  EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine forms the eigenvectors of a NONSYMMETRIC
C     TRIDIAGONAL matrix by back transforming those of the
C     corresponding symmetric matrix determined by  FIGI.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, T and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix T.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        T contains the nonsymmetric matrix.  Its subdiagonal is
C          stored in the last N-1 positions of the first column,
C          its diagonal in the N positions of the second column,
C          and its superdiagonal in the first N-1 positions of
C          the third column.  T(1,1) and T(N,3) are arbitrary.
C          T is a two-dimensional REAL array, dimensioned T(NM,3).
C
C        E contains the subdiagonal elements of the symmetric
C          matrix in its last N-1 positions.  E(1) is arbitrary.
C          E is a one-dimensional REAL array, dimensioned E(N).
C
C        M is the number of eigenvectors to be back transformed.
C          M is an INTEGER variable.
C
C        Z contains the eigenvectors to be back transformed
C          in its first M columns.  Z is a two-dimensional REAL
C          array, dimensioned Z(NM,M).
C
C     On OUTPUT
C
C        T is unaltered.
C
C        E is destroyed.
C
C        Z contains the transformed eigenvectors in its first M columns.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          2*N+I      if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
C                     In this case, the symmetric matrix is not similar
C                     to the original matrix, and the eigenvectors
C                     cannot be found by this program.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BAKVEC
C
      INTEGER I,J,M,N,NM,IERR
      REAL T(NM,3),E(*),Z(NM,*)
C
C***FIRST EXECUTABLE STATEMENT  BAKVEC
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      E(1) = 1.0E0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
         IF (E(I) .NE. 0.0E0) GO TO 80
         IF (T(I,1) .NE. 0.0E0 .OR. T(I-1,3) .NE. 0.0E0) GO TO 1000
         E(I) = 1.0E0
         GO TO 100
   80    E(I) = E(I-1) * E(I) / T(I-1,3)
  100 CONTINUE
C
      DO 120 J = 1, M
C
         DO 120 I = 2, N
         Z(I,J) = Z(I,J) * E(I)
  120 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- EIGENVECTORS CANNOT BE
C                FOUND BY THIS PROGRAM ..........
 1000 IERR = 2 * N + I
 1001 RETURN
      END
