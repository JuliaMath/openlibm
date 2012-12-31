*DECK IMTQL1
      SUBROUTINE IMTQL1 (N, D, E, IERR)
C***BEGIN PROLOGUE  IMTQL1
C***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix
C            using the implicit QL method.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (IMTQL1-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure IMTQL1,
C     NUM. MATH. 12, 377-383(1968) by Martin and Wilkinson,
C     as modified in NUM. MATH. 15, 450(1970) by Dubrulle.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     This subroutine finds the eigenvalues of a SYMMETRIC
C     TRIDIAGONAL matrix by the implicit QL method.
C
C     On INPUT
C
C        N is the order of the matrix.  N is an INTEGER variable.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional REAL array, dimensioned
C          E(N).
C
C      On OUTPUT
C
C        D contains the eigenvalues in ascending order.  If an error
C          exit is made, the eigenvalues are correct and ordered for
C          indices 1, 2, ..., IERR-1, but may not be the smallest
C          eigenvalues.
C
C        E has been destroyed.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues should be correct for indices
C                     1, 2, ..., IERR-1.  These eigenvalues are
C                     ordered, but are not necessarily the smallest.
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
C***END PROLOGUE  IMTQL1
C
      INTEGER I,J,L,M,N,II,MML,IERR
      REAL D(*),E(*)
      REAL B,C,F,G,P,R,S,S1,S2
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  IMTQL1
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = 0.0E0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            S1 = ABS(D(M)) + ABS(D(M+1))
            S2 = S1 + ABS(E(M))
            IF (S2 .EQ. S1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (D(L+1) - P) / (2.0E0 * E(L))
         R = PYTHAG(G,1.0E0)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
         S = 1.0E0
         C = 1.0E0
         P = 0.0E0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            IF (ABS(F) .LT. ABS(G)) GO TO 150
            C = G / F
            R = SQRT(C*C+1.0E0)
            E(I+1) = F * R
            S = 1.0E0 / R
            C = C * S
            GO TO 160
  150       S = F / G
            R = SQRT(S*S+1.0E0)
            E(I+1) = G * R
            C = 1.0E0 / R
            S = S * C
  160       G = D(I+1) - P
            R = (D(I) - G) * S + 2.0E0 * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = 0.0E0
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
