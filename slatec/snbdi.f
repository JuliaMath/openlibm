*DECK SNBDI
      SUBROUTINE SNBDI (ABE, LDA, N, ML, MU, IPVT, DET)
C***BEGIN PROLOGUE  SNBDI
C***PURPOSE  Compute the determinant of a band matrix using the factors
C            computed by SNBCO or SNBFA.
C***LIBRARY   SLATEC
C***CATEGORY  D3A2
C***TYPE      SINGLE PRECISION (SNBDI-S, DNBDI-D, CNBDI-C)
C***KEYWORDS  BANDED, DETERMINANT, LINEAR EQUATIONS, NONSYMMETRIC
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C     SNBDI computes the determinant of a band matrix
C     using the factors computed by SNBCO or SNBFA.
C     If the inverse is needed, use SNBSL  N  times.
C
C     On Entry
C
C        ABE     REAL(LDA, NC)
C                the output from SNBCO or SNBFA.
C                NC must be .GE. 2*ML+MU+1 .
C
C        LDA     INTEGER
C                the leading dimension of the array  ABE .
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
C                the pivot vector from SNBCO or SNBFA.
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
C   800725  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SNBDI
      INTEGER LDA,N,ML,MU,IPVT(*)
      REAL ABE(LDA,*),DET(2)
C
      REAL TEN
      INTEGER I
C***FIRST EXECUTABLE STATEMENT  SNBDI
      DET(1) = 1.0E0
      DET(2) = 0.0E0
      TEN = 10.0E0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABE(I,ML+1)*DET(1)
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
