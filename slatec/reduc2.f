*DECK REDUC2
      SUBROUTINE REDUC2 (NM, N, A, B, DL, IERR)
C***BEGIN PROLOGUE  REDUC2
C***PURPOSE  Reduce a certain generalized symmetric eigenproblem to a
C            standard symmetric eigenproblem using Cholesky
C            factorization.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1C
C***TYPE      SINGLE PRECISION (REDUC2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure REDUC2,
C     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     This subroutine reduces the generalized SYMMETRIC eigenproblems
C     ABx=(LAMBDA)x OR BAy=(LAMBDA)y, where B is POSITIVE DEFINITE,
C     to the standard symmetric eigenproblem using the Cholesky
C     factorization of B.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and B, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrices A and B.  If the Cholesky
C          factor L of B is already available, N should be prefixed
C          with a minus sign.  N is an INTEGER variable.
C
C        A and B contain the real symmetric input matrices.  Only
C          the full upper triangles of the matrices need be supplied.
C          If N is negative, the strict lower triangle of B contains,
C          instead, the strict lower triangle of its Cholesky factor L.
C          A and B are two-dimensional REAL arrays, dimensioned A(NM,N)
C          and B(NM,N).
C
C       DL contains, if N is negative, the diagonal elements of L.
C          DL is a one-dimensional REAL array, dimensioned DL(N).
C
C     On Output
C
C        A contains in its full lower triangle the full lower triangle
C          of the symmetric matrix derived from the reduction to the
C          standard form.  The strict upper triangle of A is unaltered.
C
C        B contains in its strict lower triangle the strict lower
C          triangle of its Cholesky factor L.  The full upper triangle
C          of B is unaltered.
C
C        DL contains the diagonal elements of L.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          7*N+1      if B is not positive definite.
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
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  REDUC2
C
      INTEGER I,J,K,N,I1,J1,NM,NN,IERR
      REAL A(NM,*),B(NM,*),DL(*)
      REAL X,Y
C
C***FIRST EXECUTABLE STATEMENT  REDUC2
      IERR = 0
      NN = ABS(N)
      IF (N .LT. 0) GO TO 100
C     .......... FORM L IN THE ARRAYS B AND DL ..........
      DO 80 I = 1, N
         I1 = I - 1
C
         DO 80 J = I, N
            X = B(I,J)
            IF (I .EQ. 1) GO TO 40
C
            DO 20 K = 1, I1
   20       X = X - B(I,K) * B(J,K)
C
   40       IF (J .NE. I) GO TO 60
            IF (X .LE. 0.0E0) GO TO 1000
            Y = SQRT(X)
            DL(I) = Y
            GO TO 80
   60       B(J,I) = X / Y
   80 CONTINUE
C     .......... FORM THE LOWER TRIANGLE OF A*L
C                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  100 DO 200 I = 1, NN
         I1 = I + 1
C
         DO 200 J = 1, I
            X = A(J,I) * DL(J)
            IF (J .EQ. I) GO TO 140
            J1 = J + 1
C
            DO 120 K = J1, I
  120       X = X + A(K,I) * B(K,J)
C
  140       IF (I .EQ. NN) GO TO 180
C
            DO 160 K = I1, NN
  160       X = X + A(I,K) * B(K,J)
C
  180       A(I,J) = X
  200 CONTINUE
C     .......... PRE-MULTIPLY BY TRANSPOSE(L) AND OVERWRITE ..........
      DO 300 I = 1, NN
         I1 = I + 1
         Y = DL(I)
C
         DO 300 J = 1, I
            X = Y * A(I,J)
            IF (I .EQ. NN) GO TO 280
C
            DO 260 K = I1, NN
  260       X = X + A(K,J) * B(K,I)
C
  280       A(I,J) = X
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
 1000 IERR = 7 * N + 1
 1001 RETURN
      END
