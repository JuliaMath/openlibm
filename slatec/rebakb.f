*DECK REBAKB
      SUBROUTINE REBAKB (NM, N, B, DL, M, Z)
C***BEGIN PROLOGUE  REBAKB
C***PURPOSE  Form the eigenvectors of a generalized symmetric
C            eigensystem from the eigenvectors of derived matrix output
C            from REDUC2.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (REBAKB-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure REBAKB,
C     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
C
C     This subroutine forms the eigenvectors of a generalized
C     SYMMETRIC eigensystem by back transforming those of the
C     derived symmetric matrix determined by  REDUC2.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, B and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix system.  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        B contains information about the similarity transformation
C          (Cholesky decomposition) used in the reduction by  REDUC2
C          in its strict lower triangle.  B is a two-dimensional
C          REAL array, dimensioned B(NM,N).
C
C        DL contains further information about the transformation.
C          DL is a one-dimensional REAL array, dimensioned DL(N).
C
C        M is the number of eigenvectors to be back transformed.
C          M is an INTEGER variable.
C
C        Z contains the eigenvectors to be back transformed in its
C          first M columns.  Z is a two-dimensional REAL array
C          dimensioned Z(NM,M).
C
C     On Output
C
C        Z contains the transformed eigenvectors in its first
C          M columns.
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
C***END PROLOGUE  REBAKB
C
      INTEGER I,J,K,M,N,I1,II,NM
      REAL B(NM,*),DL(*),Z(NM,*)
      REAL X
C
C***FIRST EXECUTABLE STATEMENT  REBAKB
      IF (M .EQ. 0) GO TO 200
C
      DO 100 J = 1, M
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
         DO 100 II = 1, N
            I1 = N - II
            I = I1 + 1
            X = DL(I) * Z(I,J)
            IF (I .EQ. 1) GO TO 80
C
            DO 60 K = 1, I1
   60       X = X + B(I,K) * Z(K,J)
C
   80       Z(I,J) = X
  100 CONTINUE
C
  200 RETURN
      END
