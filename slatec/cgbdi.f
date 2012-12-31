*DECK CGBDI
      SUBROUTINE CGBDI (ABD, LDA, N, ML, MU, IPVT, DET)
C***BEGIN PROLOGUE  CGBDI
C***PURPOSE  Compute the determinant of a complex band matrix using the
C            factors from CGBCO or CGBFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D3C2
C***TYPE      COMPLEX (SGBDI-S, DGBDI-D, CGBDI-C)
C***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
C             MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGBDI computes the determinant of a band matrix
C     using the factors computed by CGBCO or CGBFA.
C     If the inverse is needed, use CGBSL  N  times.
C
C     On Entry
C
C        ABD     COMPLEX(LDA, N)
C                the output from CGBCO or CGBFA.
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
C                the pivot vector from CGBCO or CGBFA.
C
C     On Return
C
C        DET     COMPLEX(2)
C                determinant of original matrix.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0
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
C***END PROLOGUE  CGBDI
      INTEGER LDA,N,ML,MU,IPVT(*)
      COMPLEX ABD(LDA,*),DET(2)
C
      REAL TEN
      INTEGER I,M
      COMPLEX ZDUM
      REAL CABS1
C
      CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CGBDI
      M = ML + MU + 1
      DET(1) = (1.0E0,0.0E0)
      DET(2) = (0.0E0,0.0E0)
      TEN = 10.0E0
      DO 50 I = 1, N
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)
         DET(1) = ABD(M,I)*DET(1)
         IF (CABS1(DET(1)) .EQ. 0.0E0) GO TO 60
   10    IF (CABS1(DET(1)) .GE. 1.0E0) GO TO 20
            DET(1) = CMPLX(TEN,0.0E0)*DET(1)
            DET(2) = DET(2) - (1.0E0,0.0E0)
         GO TO 10
   20    CONTINUE
   30    IF (CABS1(DET(1)) .LT. TEN) GO TO 40
            DET(1) = DET(1)/CMPLX(TEN,0.0E0)
            DET(2) = DET(2) + (1.0E0,0.0E0)
         GO TO 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
