*DECK ELMBAK
      SUBROUTINE ELMBAK (NM, LOW, IGH, A, INT, M, Z)
C***BEGIN PROLOGUE  ELMBAK
C***PURPOSE  Form the eigenvectors of a real general matrix from the
C            eigenvectors of the upper Hessenberg matrix output from
C            ELMHES.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (ELMBAK-S, COMBAK-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMBAK,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a REAL GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  ELMHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  BALANC.  If  BALANC  has not been
C          used, set LOW=1 and IGH equal to the order of the matrix.
C
C        A contains the multipliers which were used in the reduction
C          by  ELMHES  in its lower triangle below the subdiagonal.
C          A is a two-dimensional REAL array, dimensioned A(NM,IGH).
C
C        INT contains information on the rows and columns interchanged
C          in the reduction by  ELMHES.  Only elements LOW through IGH
C          are used.  INT is a one-dimensional INTEGER array,
C          dimensioned INT(IGH).
C
C        M is the number of columns of Z to be back transformed.
C          M is an INTEGER variable.
C
C        Z contains the real and imaginary parts of the eigenvectors
C          to be back transformed in its first M columns.  Z is a
C          two-dimensional REAL array, dimensioned Z(NM,M).
C
C     On OUTPUT
C
C        Z contains the real and imaginary parts of the transformed
C          eigenvectors in its first M columns.
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
C***END PROLOGUE  ELMBAK
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,*),Z(NM,*)
      REAL X
      INTEGER INT(*)
C
C***FIRST EXECUTABLE STATEMENT  ELMBAK
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         MP1 = MP + 1
C
         DO 110 I = MP1, IGH
            X = A(I,MP-1)
            IF (X .EQ. 0.0E0) GO TO 110
C
            DO 100 J = 1, M
  100       Z(I,J) = Z(I,J) + X * Z(MP,J)
C
  110    CONTINUE
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = 1, M
            X = Z(I,J)
            Z(I,J) = Z(MP,J)
            Z(MP,J) = X
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
