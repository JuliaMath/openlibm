*DECK DQCHEB
      SUBROUTINE DQCHEB (X, FVAL, CHEB12, CHEB24)
C***BEGIN PROLOGUE  DQCHEB
C***SUBSIDIARY
C***PURPOSE  This routine computes the CHEBYSHEV series expansion
C            of degrees 12 and 24 of a function using A
C            FAST FOURIER TRANSFORM METHOD
C            F(X) = SUM(K=1,..,13) (CHEB12(K)*T(K-1,X)),
C            F(X) = SUM(K=1,..,25) (CHEB24(K)*T(K-1,X)),
C            Where T(K,X) is the CHEBYSHEV POLYNOMIAL OF DEGREE K.
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QCHEB-S, DQCHEB-D)
C***KEYWORDS  CHEBYSHEV SERIES EXPANSION, FAST FOURIER TRANSFORM
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Chebyshev Series Expansion
C        Standard Fortran Subroutine
C        Double precision version
C
C        PARAMETERS
C          ON ENTRY
C           X      - Double precision
C                    Vector of dimension 11 containing the
C                    Values COS(K*PI/24), K = 1, ..., 11
C
C           FVAL   - Double precision
C                    Vector of dimension 25 containing the
C                    function values at the points
C                    (B+A+(B-A)*COS(K*PI/24))/2, K = 0, ...,24,
C                    where (A,B) is the approximation interval.
C                    FVAL(1) and FVAL(25) are divided by two
C                    (these values are destroyed at output).
C
C          ON RETURN
C           CHEB12 - Double precision
C                    Vector of dimension 13 containing the
C                    CHEBYSHEV coefficients for degree 12
C
C           CHEB24 - Double precision
C                    Vector of dimension 25 containing the
C                    CHEBYSHEV Coefficients for degree 24
C
C***SEE ALSO  DQC25C, DQC25F, DQC25S
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   810101  DATE WRITTEN
C   830518  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQCHEB
C
      DOUBLE PRECISION ALAM,ALAM1,ALAM2,CHEB12,CHEB24,FVAL,PART1,PART2,
     1  PART3,V,X
      INTEGER I,J
C
      DIMENSION CHEB12(13),CHEB24(25),FVAL(25),V(12),X(11)
C
C***FIRST EXECUTABLE STATEMENT  DQCHEB
      DO 10 I=1,12
        J = 26-I
        V(I) = FVAL(I)-FVAL(J)
        FVAL(I) = FVAL(I)+FVAL(J)
   10 CONTINUE
      ALAM1 = V(1)-V(9)
      ALAM2 = X(6)*(V(3)-V(7)-V(11))
      CHEB12(4) = ALAM1+ALAM2
      CHEB12(10) = ALAM1-ALAM2
      ALAM1 = V(2)-V(8)-V(10)
      ALAM2 = V(4)-V(6)-V(12)
      ALAM = X(3)*ALAM1+X(9)*ALAM2
      CHEB24(4) = CHEB12(4)+ALAM
      CHEB24(22) = CHEB12(4)-ALAM
      ALAM = X(9)*ALAM1-X(3)*ALAM2
      CHEB24(10) = CHEB12(10)+ALAM
      CHEB24(16) = CHEB12(10)-ALAM
      PART1 = X(4)*V(5)
      PART2 = X(8)*V(9)
      PART3 = X(6)*V(7)
      ALAM1 = V(1)+PART1+PART2
      ALAM2 = X(2)*V(3)+PART3+X(10)*V(11)
      CHEB12(2) = ALAM1+ALAM2
      CHEB12(12) = ALAM1-ALAM2
      ALAM = X(1)*V(2)+X(3)*V(4)+X(5)*V(6)+X(7)*V(8)
     1  +X(9)*V(10)+X(11)*V(12)
      CHEB24(2) = CHEB12(2)+ALAM
      CHEB24(24) = CHEB12(2)-ALAM
      ALAM = X(11)*V(2)-X(9)*V(4)+X(7)*V(6)-X(5)*V(8)
     1  +X(3)*V(10)-X(1)*V(12)
      CHEB24(12) = CHEB12(12)+ALAM
      CHEB24(14) = CHEB12(12)-ALAM
      ALAM1 = V(1)-PART1+PART2
      ALAM2 = X(10)*V(3)-PART3+X(2)*V(11)
      CHEB12(6) = ALAM1+ALAM2
      CHEB12(8) = ALAM1-ALAM2
      ALAM = X(5)*V(2)-X(9)*V(4)-X(1)*V(6)
     1  -X(11)*V(8)+X(3)*V(10)+X(7)*V(12)
      CHEB24(6) = CHEB12(6)+ALAM
      CHEB24(20) = CHEB12(6)-ALAM
      ALAM = X(7)*V(2)-X(3)*V(4)-X(11)*V(6)+X(1)*V(8)
     1  -X(9)*V(10)-X(5)*V(12)
      CHEB24(8) = CHEB12(8)+ALAM
      CHEB24(18) = CHEB12(8)-ALAM
      DO 20 I=1,6
        J = 14-I
        V(I) = FVAL(I)-FVAL(J)
        FVAL(I) = FVAL(I)+FVAL(J)
   20 CONTINUE
      ALAM1 = V(1)+X(8)*V(5)
      ALAM2 = X(4)*V(3)
      CHEB12(3) = ALAM1+ALAM2
      CHEB12(11) = ALAM1-ALAM2
      CHEB12(7) = V(1)-V(5)
      ALAM = X(2)*V(2)+X(6)*V(4)+X(10)*V(6)
      CHEB24(3) = CHEB12(3)+ALAM
      CHEB24(23) = CHEB12(3)-ALAM
      ALAM = X(6)*(V(2)-V(4)-V(6))
      CHEB24(7) = CHEB12(7)+ALAM
      CHEB24(19) = CHEB12(7)-ALAM
      ALAM = X(10)*V(2)-X(6)*V(4)+X(2)*V(6)
      CHEB24(11) = CHEB12(11)+ALAM
      CHEB24(15) = CHEB12(11)-ALAM
      DO 30 I=1,3
        J = 8-I
        V(I) = FVAL(I)-FVAL(J)
        FVAL(I) = FVAL(I)+FVAL(J)
   30 CONTINUE
      CHEB12(5) = V(1)+X(8)*V(3)
      CHEB12(9) = FVAL(1)-X(8)*FVAL(3)
      ALAM = X(4)*V(2)
      CHEB24(5) = CHEB12(5)+ALAM
      CHEB24(21) = CHEB12(5)-ALAM
      ALAM = X(8)*FVAL(2)-FVAL(4)
      CHEB24(9) = CHEB12(9)+ALAM
      CHEB24(17) = CHEB12(9)-ALAM
      CHEB12(1) = FVAL(1)+FVAL(3)
      ALAM = FVAL(2)+FVAL(4)
      CHEB24(1) = CHEB12(1)+ALAM
      CHEB24(25) = CHEB12(1)-ALAM
      CHEB12(13) = V(1)-V(3)
      CHEB24(13) = CHEB12(13)
      ALAM = 0.1D+01/0.6D+01
      DO 40 I=2,12
        CHEB12(I) = CHEB12(I)*ALAM
   40 CONTINUE
      ALAM = 0.5D+00*ALAM
      CHEB12(1) = CHEB12(1)*ALAM
      CHEB12(13) = CHEB12(13)*ALAM
      DO 50 I=2,24
        CHEB24(I) = CHEB24(I)*ALAM
   50 CONTINUE
      CHEB24(1) = 0.5D+00*ALAM*CHEB24(1)
      CHEB24(25) = 0.5D+00*ALAM*CHEB24(25)
      RETURN
      END
