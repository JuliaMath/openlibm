*DECK DGTSL
      SUBROUTINE DGTSL (N, C, D, E, B, INFO)
C***BEGIN PROLOGUE  DGTSL
C***PURPOSE  Solve a tridiagonal linear system.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A2A
C***TYPE      DOUBLE PRECISION (SGTSL-S, DGTSL-D, CGTSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL
C***AUTHOR  Dongarra, J., (ANL)
C***DESCRIPTION
C
C     DGTSL given a general tridiagonal matrix and a right hand
C     side will find the solution.
C
C     On Entry
C
C        N       INTEGER
C                is the order of the tridiagonal matrix.
C
C        C       DOUBLE PRECISION(N)
C                is the subdiagonal of the tridiagonal matrix.
C                C(2) through C(N) should contain the subdiagonal.
C                On output C is destroyed.
C
C        D       DOUBLE PRECISION(N)
C                is the diagonal of the tridiagonal matrix.
C                On output D is destroyed.
C
C        E       DOUBLE PRECISION(N)
C                is the superdiagonal of the tridiagonal matrix.
C                E(1) through E(N-1) should contain the superdiagonal.
C                On output E is destroyed.
C
C        B       DOUBLE PRECISION(N)
C                is the right hand side vector.
C
C     On Return
C
C        B       is the solution vector.
C
C        INFO    INTEGER
C                = 0 normal value.
C                = K if the K-th element of the diagonal becomes
C                    exactly zero.  The subroutine returns when
C                    this is detected.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGTSL
      INTEGER N,INFO
      DOUBLE PRECISION C(*),D(*),E(*),B(*)
C
      INTEGER K,KB,KP1,NM1,NM2
      DOUBLE PRECISION T
C***FIRST EXECUTABLE STATEMENT  DGTSL
         INFO = 0
         C(1) = D(1)
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 40
            D(1) = E(1)
            E(1) = 0.0D0
            E(N) = 0.0D0
C
            DO 30 K = 1, NM1
               KP1 = K + 1
C
C              FIND THE LARGEST OF THE TWO ROWS
C
               IF (ABS(C(KP1)) .LT. ABS(C(K))) GO TO 10
C
C                 INTERCHANGE ROW
C
                  T = C(KP1)
                  C(KP1) = C(K)
                  C(K) = T
                  T = D(KP1)
                  D(KP1) = D(K)
                  D(K) = T
                  T = E(KP1)
                  E(KP1) = E(K)
                  E(K) = T
                  T = B(KP1)
                  B(KP1) = B(K)
                  B(K) = T
   10          CONTINUE
C
C              ZERO ELEMENTS
C
               IF (C(K) .NE. 0.0D0) GO TO 20
                  INFO = K
                  GO TO 100
   20          CONTINUE
               T = -C(KP1)/C(K)
               C(KP1) = D(KP1) + T*D(K)
               D(KP1) = E(KP1) + T*E(K)
               E(KP1) = 0.0D0
               B(KP1) = B(KP1) + T*B(K)
   30       CONTINUE
   40    CONTINUE
         IF (C(N) .NE. 0.0D0) GO TO 50
            INFO = N
         GO TO 90
   50    CONTINUE
C
C           BACK SOLVE
C
            NM2 = N - 2
            B(N) = B(N)/C(N)
            IF (N .EQ. 1) GO TO 80
               B(NM1) = (B(NM1) - D(NM1)*B(N))/C(NM1)
               IF (NM2 .LT. 1) GO TO 70
               DO 60 KB = 1, NM2
                  K = NM2 - KB + 1
                  B(K) = (B(K) - D(K)*B(K+1) - E(K)*B(K+2))/C(K)
   60          CONTINUE
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
C
      RETURN
      END
