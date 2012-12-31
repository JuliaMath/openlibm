*DECK HTRIBK
      SUBROUTINE HTRIBK (NM, N, AR, AI, TAU, M, ZR, ZI)
C***BEGIN PROLOGUE  HTRIBK
C***PURPOSE  Form the eigenvectors of a complex Hermitian matrix from
C            the eigenvectors of a real symmetric tridiagonal matrix
C            output from HTRIDI.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (HTRIBK-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of a complex analogue of
C     the ALGOL procedure TRBAK1, NUM. MATH. 11, 181-195(1968)
C     by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine forms the eigenvectors of a COMPLEX HERMITIAN
C     matrix by back transforming those of the corresponding
C     real symmetric tridiagonal matrix determined by  HTRIDI.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, AR, AI, ZR, and ZI, as declared in the
C          calling program dimension statement.  NM is an INTEGER
C          variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        AR and AI contain some information about the unitary
C          transformations used in the reduction by  HTRIDI  in the
C          strict lower triangle of AR and the full lower triangle of
C          AI.  The remaining upper parts of the matrices are arbitrary.
C          AR and AI are two-dimensional REAL arrays, dimensioned
C          AR(NM,N) and AI(NM,N).
C
C        TAU contains further information about the transformations.
C          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
C
C        M is the number of eigenvectors to be back transformed.
C          M is an INTEGER variable.
C
C       ZR contains the eigenvectors to be back transformed in its first
C          M columns.  The contents of ZI are immaterial.  ZR and ZI are
C          two-dimensional REAL arrays, dimensioned ZR(NM,M) and
C          ZI(NM,M).
C
C     On OUTPUT
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the transformed eigenvectors in their first M columns.
C
C     Note that the last component of each returned vector
C     is real and that vector Euclidean norms are preserved.
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
C***END PROLOGUE  HTRIBK
C
      INTEGER I,J,K,L,M,N,NM
      REAL AR(NM,*),AI(NM,*),TAU(2,*),ZR(NM,*),ZI(NM,*)
      REAL H,S,SI
C
C***FIRST EXECUTABLE STATEMENT  HTRIBK
      IF (M .EQ. 0) GO TO 200
C     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C                TRIDIAGONAL MATRIX. ..........
      DO 50 K = 1, N
C
         DO 50 J = 1, M
            ZI(K,J) = -ZR(K,J) * TAU(2,K)
            ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
C
      IF (N .EQ. 1) GO TO 200
C     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. 0.0E0) GO TO 140
C
         DO 130 J = 1, M
            S = 0.0E0
            SI = 0.0E0
C
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
  110       CONTINUE
C     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
            S = (S / H) / H
            SI = (SI / H) / H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
