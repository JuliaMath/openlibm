*DECK SGBDI
      SUBROUTINE SGBDI (ABD, LDA, N, ML, MU, IPVT, DET)
C***BEGIN PROLOGUE  SGBDI
C***PURPOSE  Compute the determinant of a band matrix using the factors
C            computed by SGBCO or SGBFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D3A2
C***TYPE      SINGLE PRECISION (SGBDI-S, DGBDI-D, CGBDI-C)
C***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
C             MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     SGBDI computes the determinant of a band matrix
C     using the factors computed by SBGCO or SGBFA.
C     If the inverse is needed, use SGBSL  N  times.
C
C     On Entry
C
C        ABD     REAL(LDA, N)
C                the output from SBGCO or SGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from SBGCO or SGBFA.
C
C     On Return
C
C        DET     REAL(2)
C                determinant of original matrix.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                or  DET(1) = 0.0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SGBDI
      INTEGER LDA,N,ML,MU,IPVT(*)
      REAL ABD(LDA,*),DET(2)
C
      REAL TEN
      INTEGER I,M
C***FIRST EXECUTABLE STATEMENT  SGBDI
      M = ML + MU + 1
      DET(1) = 1.0E0
      DET(2) = 0.0E0
      TEN = 10.0E0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABD(M,I)*DET(1)
         IF (DET(1) .EQ. 0.0E0) GO TO 60
   10    IF (ABS(DET(1)) .GE. 1.0E0) GO TO 20
            DET(1) = TEN*DET(1)
            DET(2) = DET(2) - 1.0E0
         GO TO 10
   20    CONTINUE
   30    IF (ABS(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/TEN
            DET(2) = DET(2) + 1.0E0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
