*DECK DPTSL
      SUBROUTINE DPTSL (N, D, E, B)
C***BEGIN PROLOGUE  DPTSL
C***PURPOSE  Solve a positive definite tridiagonal linear system.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2B2A
C***TYPE      DOUBLE PRECISION (SPTSL-S, DPTSL-D, CPTSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE,
C             TRIDIAGONAL
C***AUTHOR  Dongarra, J., (ANL)
C***DESCRIPTION
C
C     DPTSL, given a positive definite symmetric tridiagonal matrix and
C     a right hand side, will find the solution.
C
C     On Entry
C
C        N        INTEGER
C                 is the order of the tridiagonal matrix.
C
C        D        DOUBLE PRECISION(N)
C                 is the diagonal of the tridiagonal matrix.
C                 On output D is destroyed.
C
C        E        DOUBLE PRECISION(N)
C                 is the offdiagonal of the tridiagonal matrix.
C                 E(1) through E(N-1) should contain the
C                 offdiagonal.
C
C        B        DOUBLE PRECISION(N)
C                 is the right hand side vector.
C
C     On Return
C
C        B        contains the solution.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890505  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPTSL
      INTEGER N
      DOUBLE PRECISION D(*),E(*),B(*)
C
      INTEGER K,KBM1,KE,KF,KP1,NM1,NM1D2
      DOUBLE PRECISION T1,T2
C
C     CHECK FOR 1 X 1 CASE
C
C***FIRST EXECUTABLE STATEMENT  DPTSL
      IF (N .NE. 1) GO TO 10
         B(1) = B(1)/D(1)
      GO TO 70
   10 CONTINUE
         NM1 = N - 1
         NM1D2 = NM1/2
         IF (N .EQ. 2) GO TO 30
            KBM1 = N - 1
C
C           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF
C           SUPERDIAGONAL
C
            DO 20 K = 1, NM1D2
               T1 = E(K)/D(K)
               D(K+1) = D(K+1) - T1*E(K)
               B(K+1) = B(K+1) - T1*B(K)
               T2 = E(KBM1)/D(KBM1+1)
               D(KBM1) = D(KBM1) - T2*E(KBM1)
               B(KBM1) = B(KBM1) - T2*B(KBM1+1)
               KBM1 = KBM1 - 1
   20       CONTINUE
   30    CONTINUE
         KP1 = NM1D2 + 1
C
C        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER
C
         IF (MOD(N,2) .NE. 0) GO TO 40
            T1 = E(KP1)/D(KP1)
            D(KP1+1) = D(KP1+1) - T1*E(KP1)
            B(KP1+1) = B(KP1+1) - T1*B(KP1)
            KP1 = KP1 + 1
   40    CONTINUE
C
C        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP
C        AND BOTTOM
C
         B(KP1) = B(KP1)/D(KP1)
         IF (N .EQ. 2) GO TO 60
            K = KP1 - 1
            KE = KP1 + NM1D2 - 1
            DO 50 KF = KP1, KE
               B(K) = (B(K) - E(K)*B(K+1))/D(K)
               B(KF+1) = (B(KF+1) - E(KF)*B(KF))/D(KF+1)
               K = K - 1
   50       CONTINUE
   60    CONTINUE
         IF (MOD(N,2) .EQ. 0) B(1) = (B(1) - E(1)*B(2))/D(1)
   70 CONTINUE
      RETURN
      END
