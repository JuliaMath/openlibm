*DECK ORTRAN
      SUBROUTINE ORTRAN (NM, N, LOW, IGH, A, ORT, Z)
C***BEGIN PROLOGUE  ORTRAN
C***PURPOSE  Accumulate orthogonal similarity transformations in the
C            reduction of real general matrix by ORTHES.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (ORTRAN-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ORTRANS,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine accumulates the orthogonal similarity
C     transformations used in the reduction of a REAL GENERAL
C     matrix to upper Hessenberg form by  ORTHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  BALANC.  If  BALANC  has not been
C          used, set LOW=1 and IGH equal to the order of the matrix, N.
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
C     On OUTPUT
C
C        Z contains the transformation matrix produced in the reduction
C          by  ORTHES  to the upper Hessenberg form.  Z is a two-
C          dimensional REAL array, dimensioned Z(NM,N).
C
C        ORT has been used for temporary storage as is not restored.
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
C***END PROLOGUE  ORTRAN
C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,*),ORT(*),Z(NM,*)
      REAL G
C
C     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
C***FIRST EXECUTABLE STATEMENT  ORTRAN
      DO 80 I = 1, N
C
         DO 60 J = 1, N
   60    Z(I,J) = 0.0E0
C
         Z(I,I) = 1.0E0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0E0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = MP, IGH
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
