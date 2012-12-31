*DECK HTRID3
      SUBROUTINE HTRID3 (NM, N, A, D, E, E2, TAU)
C***BEGIN PROLOGUE  HTRID3
C***PURPOSE  Reduce a complex Hermitian (packed) matrix to a real
C            symmetric tridiagonal matrix by unitary similarity
C            transformations.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B1
C***TYPE      SINGLE PRECISION (HTRID3-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of a complex analogue of
C     the ALGOL procedure TRED3, NUM. MATH. 11, 181-195(1968)
C     by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a COMPLEX HERMITIAN matrix, stored as
C     a single square array, to a real symmetric tridiagonal matrix
C     using unitary similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, A, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        A contains the lower triangle of the complex Hermitian input
C          matrix.  The real parts of the matrix elements are stored
C          in the full lower triangle of A, and the imaginary parts
C          are stored in the transposed positions of the strict upper
C          triangle of A.  No storage is required for the zero
C          imaginary parts of the diagonal elements.  A is a two-
C          dimensional REAL array, dimensioned A(NM,N).
C
C     On OUTPUT
C
C        A contains some information about the unitary transformations
C          used in the reduction.
C
C        D contains the diagonal elements of the real symmetric
C          tridiagonal matrix.  D is a one-dimensional REAL array,
C          dimensioned D(N).
C
C        E contains the subdiagonal elements of the real tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C          E is a one-dimensional REAL array, dimensioned E(N).
C
C        E2 contains the squares of the corresponding elements of E.
C          E2(1) is set to zero.  E2 may coincide with E if the squares
C          are not needed.  E2 is a one-dimensional REAL array,
C          dimensioned E2(N).
C
C        TAU contains further information about the transformations.
C          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  HTRID3
C
      INTEGER I,J,K,L,N,II,NM,JM1,JP1
      REAL A(NM,*),D(*),E(*),E2(*),TAU(2,*)
      REAL F,G,H,FI,GI,HH,SI,SCALE
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  HTRID3
      TAU(1,N) = 1.0E0
      TAU(2,N) = 0.0E0
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K)) + ABS(A(K,I))
C
         IF (SCALE .NE. 0.0E0) GO TO 140
         TAU(1,L) = 1.0E0
         TAU(2,L) = 0.0E0
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            A(K,I) = A(K,I) / SCALE
            H = H + A(I,K) * A(I,K) + A(K,I) * A(K,I)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
         F = PYTHAG(A(I,L),A(L,I))
C     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
         IF (F .EQ. 0.0E0) GO TO 160
         TAU(1,L) = (A(L,I) * TAU(2,I) - A(I,L) * TAU(1,I)) / F
         SI = (A(I,L) * TAU(2,I) + A(L,I) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0E0 + G / F
         A(I,L) = G * A(I,L)
         A(L,I) = G * A(L,I)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         A(I,L) = G
  170    F = 0.0E0
C
         DO 240 J = 1, L
            G = 0.0E0
            GI = 0.0E0
            IF (J .EQ. 1) GO TO 190
            JM1 = J - 1
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, JM1
               G = G + A(J,K) * A(I,K) + A(K,J) * A(K,I)
               GI = GI - A(J,K) * A(K,I) + A(K,J) * A(I,K)
  180       CONTINUE
C
  190       G = G + A(J,J) * A(I,J)
            GI = GI - A(J,J) * A(J,I)
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + A(K,J) * A(I,K) - A(J,K) * A(K,I)
               GI = GI - A(K,J) * A(K,I) - A(J,K) * A(I,K)
  200       CONTINUE
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * A(I,J) - TAU(2,J) * A(J,I)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -A(J,I)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
            A(J,J) = A(J,J) - 2.0E0 * (F * G + FI * GI)
            IF (J .EQ. 1) GO TO 260
            JM1 = J - 1
C
            DO 250 K = 1, JM1
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
     1                         + FI * TAU(2,K) + GI * A(K,I)
               A(K,J) = A(K,J) - F * TAU(2,K) - G * A(K,I)
     1                         - FI * E(K) - GI * A(I,K)
  250       CONTINUE
C
  260    CONTINUE
C
  270    DO 280 K = 1, L
            A(I,K) = SCALE * A(I,K)
            A(K,I) = SCALE * A(K,I)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    D(I) = A(I,I)
         A(I,I) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
