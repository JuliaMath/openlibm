*DECK ELTRAN
      SUBROUTINE ELTRAN (NM, N, LOW, IGH, A, INT, Z)
C***BEGIN PROLOGUE  ELTRAN
C***PURPOSE  Accumulates the stabilized elementary similarity
C            transformations used in the reduction of a real general
C            matrix to upper Hessenberg form by ELMHES.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      SINGLE PRECISION (ELTRAN-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     This subroutine accumulates the stabilized elementary
C     similarity transformations used in the reduction of a
C     REAL GENERAL matrix to upper Hessenberg form by  ELMHES.
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
C        A contains the multipliers which were used in the reduction
C          by  ELMHES  in its lower triangle below the subdiagonal.
C          A is a two-dimensional REAL array, dimensioned A(NM,IGH).
C
C        INT contains information on the rows and columns interchanged
C          in the reduction by  ELMHES.  Only elements LOW through IGH
C          are used.  INT is a one-dimensional INTEGER array,
C          dimensioned INT(IGH).
C
C     On OUTPUT
C
C        Z contains the transformation matrix produced in the reduction
C          by  ELMHES.  Z is a two-dimensional REAL array, dimensioned
C          Z(NM,N).
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
C***END PROLOGUE  ELTRAN
C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,*),Z(NM,*)
      INTEGER INT(*)
C
C***FIRST EXECUTABLE STATEMENT  ELTRAN
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
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0E0
  130    CONTINUE
C
         Z(I,MP) = 1.0E0
  140 CONTINUE
C
  200 RETURN
      END
