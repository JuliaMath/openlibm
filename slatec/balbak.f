*DECK BALBAK
      SUBROUTINE BALBAK (NM, N, LOW, IGH, SCALE, M, Z)
C***BEGIN PROLOGUE  BALBAK
C***PURPOSE  Form the eigenvectors of a real general matrix from the
C            eigenvectors of matrix output from BALANC.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (BALBAK-S, CBABK2-C)
C***KEYWORDS  EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BALBAK,
C     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
C     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
C
C     This subroutine forms the eigenvectors of a REAL GENERAL
C     matrix by back transforming those of the corresponding
C     balanced matrix determined by  BALANC.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the number of components of the vectors in matrix Z.
C          N is an INTEGER variable.  N must be less than or equal
C          to NM.
C
C        LOW and IGH are INTEGER variables determined by  BALANC.
C
C        SCALE contains information determining the permutations and
C          scaling factors used by  BALANC.  SCALE is a one-dimensional
C          REAL array, dimensioned SCALE(N).
C
C        M is the number of columns of Z to be back transformed.
C          M is an INTEGER variable.
C
C        Z contains the real and imaginary parts of the eigen-
C          vectors to be back transformed in its first M columns.
C          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
C
C     On OUTPUT
C
C        Z contains the real and imaginary parts of the
C          transformed eigenvectors in its first M columns.
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
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
C***END PROLOGUE  BALBAK
C
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(*),Z(NM,*)
      REAL S
C
C***FIRST EXECUTABLE STATEMENT  BALBAK
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0E0/SCALE(I). ..........
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
C
  110 CONTINUE
C     ......... FOR I=LOW-1 STEP -1 UNTIL 1,
C               IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
