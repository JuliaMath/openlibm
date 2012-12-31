*DECK TRED3
      SUBROUTINE TRED3 (N, NV, A, D, E, E2)
C***BEGIN PROLOGUE  TRED3
C***PURPOSE  Reduce a real symmetric matrix stored in packed form to
C            symmetric tridiagonal matrix using orthogonal
C            transformations.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B1
C***TYPE      SINGLE PRECISION (TRED3-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED3,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a REAL SYMMETRIC matrix, stored as
C     a one-dimensional array, to a symmetric tridiagonal matrix
C     using orthogonal similarity transformations.
C
C     On Input
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C
C        NV is an INTEGER variable set equal to the dimension of the
C          array A as specified in the calling program.  NV must not
C          be less than  N*(N+1)/2.
C
C        A contains the lower triangle, stored row-wise, of the real
C          symmetric packed matrix.  A is a one-dimensional REAL
C          array, dimensioned A(NV).
C
C     On Output
C
C        A contains information about the orthogonal transformations
C          used in the reduction in its first N*(N+1)/2 positions.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is set
C          to zero.  E is a one-dimensional REAL array, dimensioned
C          E(N).
C
C        E2 contains the squares of the corresponding elements of E.
C          E2 may coincide with E if the squares are not needed.
C          E2 is a one-dimensional REAL array, dimensioned E2(N).
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
C***END PROLOGUE  TRED3
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV
      REAL A(*),D(*),E(*),E2(*)
      REAL F,G,H,HH,SCALE
C
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
C***FIRST EXECUTABLE STATEMENT  TRED3
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I * L) / 2
         H = 0.0E0
         SCALE = 0.0E0
         IF (L .LT. 1) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
  120    CONTINUE
C
         IF (SCALE .NE. 0.0E0) GO TO 140
  130    E(I) = 0.0E0
         E2(I) = 0.0E0
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 290
         F = 0.0E0
C
         DO 240 J = 1, L
            G = 0.0E0
            JK = (J * (J-1)) / 2
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, L
               JK = JK + 1
               IF (K .GT. J) JK = JK + K - 2
               G = G + A(JK) * D(K)
  180       CONTINUE
C     .......... FORM ELEMENT OF P ..........
            E(J) = G / H
            F = F + E(J) * D(J)
  240    CONTINUE
C
         HH = F / (H + H)
         JK = 0
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = D(J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F * E(K) - G * D(K)
  260    CONTINUE
C
  290    D(I) = A(IZ+1)
         A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
