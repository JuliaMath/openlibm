*DECK DNBDI
      SUBROUTINE DNBDI (ABE, LDA, N, ML, MU, IPVT, DET)
C***BEGIN PROLOGUE  DNBDI
C***PURPOSE  Compute the determinant of a band matrix using the factors
C            computed by DNBCO or DNBFA.
C***LIBRARY   SLATEC
C***CATEGORY  D3A2
C***TYPE      DOUBLE PRECISION (SNBDI-S, DNBDI-D, CNBDI-C)
C***KEYWORDS  BANDED, DETERMINANT, LINEAR EQUATIONS, NONSYMMETRIC
C***AUTHOR  Voorhees, E. A., (LANL)
C***DESCRIPTION
C
C     DNBDI computes the determinant of a band matrix
C     using the factors computed by DNBCO or DNBFA.
C     If the inverse is needed, use DNBSL  N  times.
C
C     On Entry
C
C        ABE     DOUBLE PRECISION(LDA, NC)
C                the output from DNBCO or DNBFA.
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
C                the pivot vector from DNBCO or DNBFA.
C
C     On Return
C
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                or  DET(1) = 0.0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800728  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNBDI
      INTEGER LDA,N,ML,MU,IPVT(*)
      DOUBLE PRECISION ABE(LDA,*),DET(2)
C
      DOUBLE PRECISION TEN
      INTEGER I
C***FIRST EXECUTABLE STATEMENT  DNBDI
      DET(1) = 1.0D0
      DET(2) = 0.0D0
      TEN = 10.0D0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABE(I,ML+1)*DET(1)
         IF (DET(1) .EQ. 0.0D0) GO TO 60
   10    IF (ABS(DET(1)) .GE. 1.0D0) GO TO 20
            DET(1) = TEN*DET(1)
            DET(2) = DET(2) - 1.0D0
         GO TO 10
   20    CONTINUE
   30    IF (ABS(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/TEN
            DET(2) = DET(2) + 1.0D0
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
