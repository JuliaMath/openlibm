*DECK TRBAK3
      SUBROUTINE TRBAK3 (NM, N, NV, A, M, Z)
C***BEGIN PROLOGUE  TRBAK3
C***PURPOSE  Form the eigenvectors of a real symmetric matrix from the
C            eigenvectors of a symmetric tridiagonal matrix formed
C            by TRED3.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (TRBAK3-S)
C***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRBAK3,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine forms the eigenvectors of a REAL SYMMETRIC
C     matrix by back transforming those of the corresponding
C     symmetric tridiagonal matrix determined by  TRED3.
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
C        NV is an INTEGER variable set equal to the dimension of the
C          array A as specified in the calling program.  NV must not
C          be less than  N*(N+1)/2.
C
C        A contains information about the orthogonal transformations
C          used in the reduction by  TRED3  in its first N*(N+1)/2
C          positions.  A is a one-dimensional REAL array, dimensioned
C          A(NV).
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
C     Note that TRBAK3 preserves vector Euclidean norms.
C
C     Questions and comments should be directed to b. s. Garbow,
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
C***END PROLOGUE  TRBAK3
C
      INTEGER I,J,K,L,M,N,IK,IZ,NM,NV
      REAL A(*),Z(NM,*)
      REAL H,S
C
C***FIRST EXECUTABLE STATEMENT  TRBAK3
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IZ = (I * L) / 2
         IK = IZ + I
         H = A(IK)
         IF (H .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
            IK = IZ
C
            DO 110 K = 1, L
               IK = IK + 1
               S = S + A(IK) * Z(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            IK = IZ
C
            DO 120 K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S * A(IK)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
