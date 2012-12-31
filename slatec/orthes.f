*DECK ORTHES
      SUBROUTINE ORTHES (NM, N, LOW, IGH, A, ORT)
C***BEGIN PROLOGUE  ORTHES
C***PURPOSE  Reduce a real general matrix to upper Hessenberg form
C            using orthogonal similarity transformations.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B2
C***TYPE      SINGLE PRECISION (ORTHES-S, CORTH-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure ORTHES,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     Given a REAL GENERAL matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     LOW through IGH to upper Hessenberg form by
C     orthogonal similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, A, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  BALANC.  If  BALANC  has not been
C          used, set LOW=1 and IGH equal to the order of the matrix, N.
C
C        A contains the general matrix to be reduced to upper
C          Hessenberg form.  A is a two-dimensional REAL array,
C          dimensioned A(NM,N).
C
C     On OUTPUT
C
C        A contains the upper Hessenberg matrix.  Some information about
C          the orthogonal transformations used in the reduction
C          is stored in the remaining triangle under the Hessenberg
C          matrix.
C
C        ORT contains further information about the orthogonal trans-
C          formations used in the reduction.  Only elements LOW+1
C          through IGH are used.  ORT is a one-dimensional REAL array,
C          dimensioned ORT(IGH).
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
C***END PROLOGUE  ORTHES
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL A(NM,*),ORT(*)
      REAL F,G,H,SCALE
C
C***FIRST EXECUTABLE STATEMENT  ORTHES
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0E0
         ORT(M) = 0.0E0
         SCALE = 0.0E0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(A(I,M-1))
C
         IF (SCALE .EQ. 0.0E0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
C
         G = -SIGN(SQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
C     .......... FORM (I-(U*UT)/H) * A ..........
         DO 130 J = M, N
            F = 0.0E0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
C
            F = F / H
C
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            F = 0.0E0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
C
            F = F / H
C
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
C
  160    CONTINUE
C
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
C
  200 RETURN
      END
