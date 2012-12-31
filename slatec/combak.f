*DECK COMBAK
      SUBROUTINE COMBAK (NM, LOW, IGH, AR, AI, INT, M, ZR, ZI)
C***BEGIN PROLOGUE  COMBAK
C***PURPOSE  Form the eigenvectors of a complex general matrix from the
C            eigenvectors of a upper Hessenberg matrix output from
C            COMHES.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      COMPLEX (ELMBAK-S, COMBAK-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure COMBAK,
C     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a COMPLEX GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  COMHES.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, AR, AI, ZR and ZI, as declared in the
C          calling program dimension statement.  NM is an INTEGER
C          variable.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix.
C
C        AR and AI contain the multipliers which were used in the
C           reduction by  COMHES  in their lower triangles below
C           the subdiagonal.  AR and AI are two-dimensional REAL
C           arrays, dimensioned AR(NM,IGH) and AI(NM,IGH).
C
C        INT contains information on the rows and columns
C          interchanged in the reduction by  COMHES.  Only
C          elements LOW through IGH are used.  INT is a
C          one-dimensional INTEGER array, dimensioned INT(IGH).
C
C        M is the number of eigenvectors to be back transformed.
C          M is an INTEGER variable.
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the eigenvectors to be back transformed in their first M
C          columns.  ZR and ZI are two-dimensional REAL arrays,
C          dimensioned ZR(NM,M) and ZI(NM,M).
C
C     On OUTPUT
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the transformed eigenvectors in their first M columns.
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
C***END PROLOGUE  COMBAK
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,*),AI(NM,*),ZR(NM,*),ZI(NM,*)
      REAL XR,XI
      INTEGER INT(*)
C
C***FIRST EXECUTABLE STATEMENT  COMBAK
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
            XR = AR(I,MP-1)
            XI = AI(I,MP-1)
            IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 110
C
            DO 100 J = 1, M
               ZR(I,J) = ZR(I,J) + XR * ZR(MP,J) - XI * ZI(MP,J)
               ZI(I,J) = ZI(I,J) + XR * ZI(MP,J) + XI * ZR(MP,J)
  100       CONTINUE
C
  110    CONTINUE
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = 1, M
            XR = ZR(I,J)
            ZR(I,J) = ZR(MP,J)
            ZR(MP,J) = XR
            XI = ZI(I,J)
            ZI(I,J) = ZI(MP,J)
            ZI(MP,J) = XI
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
