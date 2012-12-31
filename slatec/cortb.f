*DECK CORTB
      SUBROUTINE CORTB (NM, LOW, IGH, AR, AI, ORTR, ORTI, M, ZR, ZI)
C***BEGIN PROLOGUE  CORTB
C***PURPOSE  Form the eigenvectors of a complex general matrix from
C            eigenvectors of upper Hessenberg matrix output from
C            CORTH.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C4
C***TYPE      COMPLEX (ORTBAK-S, CORTB-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of a complex analogue of
C     the ALGOL procedure ORTBAK, NUM. MATH. 12, 349-368(1968)
C     by Martin and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     This subroutine forms the eigenvectors of a COMPLEX GENERAL
C     matrix by back transforming those of the corresponding
C     upper Hessenberg matrix determined by  CORTH.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, AR, AI, ZR, and ZI, as declared in the
C          calling program dimension statement.  NM is an INTEGER
C          variable.
C
C        LOW and IGH are two INTEGER variables determined by the
C          balancing subroutine  CBAL.  If  CBAL  has not been used,
C          set LOW=1 and IGH equal to the order of the matrix.
C
C        AR and AI contain information about the unitary trans-
C          formations used in the reduction by  CORTH  in their
C          strict lower triangles.  AR and AI are two-dimensional
C          REAL arrays, dimensioned AR(NM,IGH) and AI(NM,IGH).
C
C        ORTR and ORTI contain further information about the unitary
C          transformations used in the reduction by  CORTH.  Only
C          elements LOW through IGH are used.  ORTR and ORTI are
C          one-dimensional REAL arrays, dimensioned ORTR(IGH) and
C          ORTI(IGH).
C
C        M is the number of columns of Z=(ZR,ZI) to be back transformed.
C          M is an INTEGER variable.
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the eigenvectors to be back transformed in their first
C          M columns.  ZR and ZI are two-dimensional REAL arrays,
C          dimensioned ZR(NM,M) and ZI(NM,M).
C
C     On OUTPUT
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the transformed eigenvectors in their first M columns.
C
C        ORTR and ORTI have been altered.
C
C     Note that CORTB preserves vector Euclidean norms.
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
C***END PROLOGUE  CORTB
C
      INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
      REAL AR(NM,*),AI(NM,*),ORTR(*),ORTI(*)
      REAL ZR(NM,*),ZI(NM,*)
      REAL H,GI,GR
C
C***FIRST EXECUTABLE STATEMENT  CORTB
      IF (M .EQ. 0) GO TO 200
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = KP1, LA
         MP = LOW + IGH - MM
         IF (AR(MP,MP-1) .EQ. 0.0E0 .AND. AI(MP,MP-1) .EQ. 0.0E0)
     1      GO TO 140
C     .......... H BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         H = AR(MP,MP-1) * ORTR(MP) + AI(MP,MP-1) * ORTI(MP)
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
            ORTR(I) = AR(I,MP-1)
            ORTI(I) = AI(I,MP-1)
  100    CONTINUE
C
         DO 130 J = 1, M
            GR = 0.0E0
            GI = 0.0E0
C
            DO 110 I = MP, IGH
               GR = GR + ORTR(I) * ZR(I,J) + ORTI(I) * ZI(I,J)
               GI = GI + ORTR(I) * ZI(I,J) - ORTI(I) * ZR(I,J)
  110       CONTINUE
C
            GR = GR / H
            GI = GI / H
C
            DO 120 I = MP, IGH
               ZR(I,J) = ZR(I,J) + GR * ORTR(I) - GI * ORTI(I)
               ZI(I,J) = ZI(I,J) + GR * ORTI(I) + GI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
