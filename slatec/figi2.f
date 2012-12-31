*DECK FIGI2
      SUBROUTINE FIGI2 (NM, N, T, D, E, Z, IERR)
C***BEGIN PROLOGUE  FIGI2
C***PURPOSE  Transforms certain real non-symmetric tridiagonal matrix
C            to symmetric tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1C
C***TYPE      SINGLE PRECISION (FIGI2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     Given a NONSYMMETRIC TRIDIAGONAL matrix such that the products
C     of corresponding pairs of off-diagonal elements are all
C     non-negative, and zero only when both factors are zero, this
C     subroutine reduces it to a SYMMETRIC TRIDIAGONAL matrix
C     using and accumulating diagonal similarity transformations.
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
C     On OUTPUT
C
C        T is unaltered.
C
C        D contains the diagonal elements of the tridiagonal symmetric
C          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the tridiagonal
C          symmetric matrix in its last N-1 positions.  E(1) is not set.
C          E is a one-dimensional REAL array, dimensioned E(N).
C
C        Z contains the diagonal transformation matrix produced in the
C          symmetrization.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          N+I        if T(I,1)*T(I-1,3) is negative,
C          2*N+I      if T(I,1)*T(I-1,3) is zero with one factor
C                     non-zero.  In these cases, there does not exist
C                     a symmetrizing similarity transformation which
C                     is essential for the validity of the later
C                     eigenvector computation.
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
C***END PROLOGUE  FIGI2
C
      INTEGER I,J,N,NM,IERR
      REAL T(NM,3),D(*),E(*),Z(NM,*)
      REAL H
C
C***FIRST EXECUTABLE STATEMENT  FIGI2
      IERR = 0
C
      DO 100 I = 1, N
C
         DO 50 J = 1, N
   50    Z(I,J) = 0.0E0
C
         IF (I .EQ. 1) GO TO 70
         H = T(I,1) * T(I-1,3)
         IF (H) 900, 60, 80
   60    IF (T(I,1) .NE. 0.0E0 .OR. T(I-1,3) .NE. 0.0E0) GO TO 1000
         E(I) = 0.0E0
   70    Z(I,I) = 1.0E0
         GO TO 90
   80    E(I) = SQRT(H)
         Z(I,I) = Z(I-1,I-1) * E(I) / T(I-1,3)
   90    D(I) = T(I,2)
  100 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS NEGATIVE ..........
  900 IERR = N + I
      GO TO 1001
C     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
C                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO ..........
 1000 IERR = 2 * N + I
 1001 RETURN
      END
