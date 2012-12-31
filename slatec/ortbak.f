*DECK ORTBAK
      SUBROUTINE ORTBAK (NM, LOW, IGH, A, ORT, M, Z)
C***BEGIN PROLOGUE  ORTBAK
C***PURPOSE  Form the eigenvectors of a general real matrix from the
C            eigenvectors of the upper Hessenberg matrix output from
C            ORTHES.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (ORTBAK-S, CORTB-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ORTBAK,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a REAL GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  ORTHES.
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
C        A contains some information about the orthogonal trans-
C          formations used in the reduction to Hessenberg form by
C          ORTHES  in its strict lower triangle.  A is a two-dimensional
C          REAL array, dimensioned A(NM,IGH).
C
C        ORT contains further information about the orthogonal trans-
C          formations used in the reduction by  ORTHES.  Only elements
C          LOW through IGH are used.  ORT is a one-dimensional REAL
C          array, dimensioned ORT(IGH).
C
C        M is the number of columns of Z to be back transformed.
C          M is an INTEGER variable.
C
C        Z contains the real and imaginary parts of the eigenvectors to
C          be back transformed in its first M columns.  Z is a two-
C          dimensional REAL array, dimensioned Z(NM,M).
C
C     On OUTPUT
C
C        Z contains the real and imaginary parts of the transformed
C          eigenvectors in its first M columns.
C
C        ORT has been used for temporary storage as is not restored.
C
C     NOTE that ORTBAK preserves vector Euclidean norms.
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
C***END PROLOGUE  ORTBAK
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL A(NM,*),ORT(*),Z(NM,*)
      REAL G
C
C***FIRST EXECUTABLE STATEMENT  ORTBAK
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0E0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = 1, M
            G = 0.0E0
C
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            G = (G / ORT(MP)) / A(MP,MP-1)
C
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
