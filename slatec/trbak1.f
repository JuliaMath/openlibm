*DECK TRBAK1
      SUBROUTINE TRBAK1 (NM, N, A, E, M, Z)
C***BEGIN PROLOGUE  TRBAK1
C***PURPOSE  Form the eigenvectors of real symmetric matrix from
C            the eigenvectors of a symmetric tridiagonal matrix formed
C            by TRED1.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (TRBAK1-S)
C***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRBAK1,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine forms the eigenvectors of a REAL SYMMETRIC
C     matrix by back transforming those of the corresponding
C     symmetric tridiagonal matrix determined by  TRED1.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        A contains information about the orthogonal transformations
C          used in the reduction by  TRED1  in its strict lower
C          triangle.  A is a two-dimensional REAL array, dimensioned
C          A(NM,N).
C
C        E contains the subdiagonal elements of the tridiagonal matrix
C          in its last N-1 positions.  E(1) is arbitrary.  These
C          elements provide the remaining information about the
C          orthogonal transformations.  E is a one-dimensional REAL
C          array, dimensioned E(N).
C
C        M is the number of columns of Z to be back transformed.
C          M is an INTEGER variable.
C
C        Z contains the eigenvectors to be back transformed in its
C          first M columns.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,M).
C
C     On Output
C
C        Z contains the transformed eigenvectors in its first M columns.
C
C     Note that TRBAK1 preserves vector Euclidean norms.
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
C***END PROLOGUE  TRBAK1
C
      INTEGER I,J,K,L,M,N,NM
      REAL A(NM,*),E(*),Z(NM,*)
      REAL S
C
C***FIRST EXECUTABLE STATEMENT  TRBAK1
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IF (E(I) .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
C
            DO 110 K = 1, L
  110       S = S + A(I,K) * Z(K,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / A(I,L)) / E(I)
C
            DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(I,K)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
