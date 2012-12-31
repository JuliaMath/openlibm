*DECK BANDR
      SUBROUTINE BANDR (NM, N, MB, A, D, E, E2, MATZ, Z)
C***BEGIN PROLOGUE  BANDR
C***PURPOSE  Reduce a real symmetric band matrix to symmetric
C            tridiagonal matrix and, optionally, accumulate
C            orthogonal similarity transformations.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B1
C***TYPE      SINGLE PRECISION (BANDR-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BANDRD,
C     NUM. MATH. 12, 231-241(1968) by Schwarz.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 273-283(1971).
C
C     This subroutine reduces a REAL SYMMETRIC BAND matrix
C     to a symmetric tridiagonal matrix using and optionally
C     accumulating orthogonal similarity transformations.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        MB is the (half) band width of the matrix, defined as the
C          number of adjacent diagonals, including the principal
C          diagonal, required to specify the non-zero portion of the
C          lower triangle of the matrix.  MB is less than or equal
C          to N.  MB is an INTEGER variable.
C
C        A contains the lower triangle of the real symmetric band
C          matrix.  Its lowest subdiagonal is stored in the last
C          N+1-MB  positions of the first column, its next subdiagonal
C          in the last  N+2-MB  positions of the second column, further
C          subdiagonals similarly, and finally its principal diagonal
C          in the  N  positions of the last column.  Contents of storage
C          locations not part of the matrix are arbitrary.  A is a
C          two-dimensional REAL array, dimensioned A(NM,MB).
C
C        MATZ should be set to .TRUE. if the transformation matrix is
C          to be accumulated, and to .FALSE. otherwise.  MATZ is a
C          LOGICAL variable.
C
C     On OUTPUT
C
C        A has been destroyed, except for its last two columns which
C          contain a copy of the tridiagonal matrix.
C
C        D contains the diagonal elements of the tridiagonal matrix.
C          D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the tridiagonal
C          matrix in its last N-1 positions.  E(1) is set to zero.
C          E is a one-dimensional REAL array, dimensioned E(N).
C
C        E2 contains the squares of the corresponding elements of E.
C          E2 may coincide with E if the squares are not needed.
C          E2 is a one-dimensional REAL array, dimensioned E2(N).
C
C        Z contains the orthogonal transformation matrix produced in
C          the reduction if MATZ has been set to .TRUE.  Otherwise, Z
C          is not referenced.  Z is a two-dimensional REAL array,
C          dimensioned Z(NM,N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BANDR
C
      INTEGER J,K,L,N,R,I1,I2,J1,J2,KR,MB,MR,M1,NM,N2,R1,UGL,MAXL,MAXR
      REAL A(NM,*),D(*),E(*),E2(*),Z(NM,*)
      REAL G,U,B1,B2,C2,F1,F2,S2,DMIN,DMINRT
      LOGICAL MATZ
C
C***FIRST EXECUTABLE STATEMENT  BANDR
      DMIN = 2.0E0**(-64)
      DMINRT = 2.0E0**(-32)
C     .......... INITIALIZE DIAGONAL SCALING MATRIX ..........
      DO 30 J = 1, N
   30 D(J) = 1.0E0
C
      IF (.NOT. MATZ) GO TO 60
C
      DO 50 J = 1, N
C
         DO 40 K = 1, N
   40    Z(J,K) = 0.0E0
C
         Z(J,J) = 1.0E0
   50 CONTINUE
C
   60 M1 = MB - 1
      IF (M1 - 1) 900, 800, 70
   70 N2 = N - 2
C
      DO 700 K = 1, N2
         MAXR = MIN(M1,N-K)
C     .......... FOR R=MAXR STEP -1 UNTIL 2 DO -- ..........
         DO 600 R1 = 2, MAXR
            R = MAXR + 2 - R1
            KR = K + R
            MR = MB - R
            G = A(KR,MR)
            A(KR-1,1) = A(KR-1,MR+1)
            UGL = K
C
            DO 500 J = KR, N, M1
               J1 = J - 1
               J2 = J1 - 1
               IF (G .EQ. 0.0E0) GO TO 600
               B1 = A(J1,1) / G
               B2 = B1 * D(J1) / D(J)
               S2 = 1.0E0 / (1.0E0 + B1 * B2)
               IF (S2 .GE. 0.5E0 ) GO TO 450
               B1 = G / A(J1,1)
               B2 = B1 * D(J) / D(J1)
               C2 = 1.0E0 - S2
               D(J1) = C2 * D(J1)
               D(J) = C2 * D(J)
               F1 = 2.0E0 * A(J,M1)
               F2 = B1 * A(J1,MB)
               A(J,M1) = -B2 * (B1 * A(J,M1) - A(J,MB)) - F2 + A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J,MB) + F1) + A(J1,MB)
               A(J,MB) = B1 * (F2 - F1) + A(J,MB)
C
               DO 200 L = UGL, J2
                  I2 = MB - J + L
                  U = A(J1,I2+1) + B2 * A(J,I2)
                  A(J,I2) = -B1 * A(J1,I2+1) + A(J,I2)
                  A(J1,I2+1) = U
  200          CONTINUE
C
               UGL = J
               A(J1,1) = A(J1,1) + B2 * G
               IF (J .EQ. N) GO TO 350
               MAXL = MIN(M1,N-J1)
C
               DO 300 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = A(I1,I2) + B2 * A(I1,I2+1)
                  A(I1,I2+1) = -B1 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2) = U
  300          CONTINUE
C
               I1 = J + M1
               IF (I1 .GT. N) GO TO 350
               G = B2 * A(I1,1)
  350          IF (.NOT. MATZ) GO TO 500
C
               DO 400 L = 1, N
                  U = Z(L,J1) + B2 * Z(L,J)
                  Z(L,J) = -B1 * Z(L,J1) + Z(L,J)
                  Z(L,J1) = U
  400          CONTINUE
C
               GO TO 500
C
  450          U = D(J1)
               D(J1) = S2 * D(J)
               D(J) = S2 * U
               F1 = 2.0E0 * A(J,M1)
               F2 = B1 * A(J,MB)
               U = B1 * (F2 - F1) + A(J1,MB)
               A(J,M1) = B2 * (B1 * A(J,M1) - A(J1,MB)) + F2 - A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J1,MB) + F1) + A(J,MB)
               A(J,MB) = U
C
               DO 460 L = UGL, J2
                  I2 = MB - J + L
                  U = B2 * A(J1,I2+1) + A(J,I2)
                  A(J,I2) = -A(J1,I2+1) + B1 * A(J,I2)
                  A(J1,I2+1) = U
  460          CONTINUE
C
               UGL = J
               A(J1,1) = B2 * A(J1,1) + G
               IF (J .EQ. N) GO TO 480
               MAXL = MIN(M1,N-J1)
C
               DO 470 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = B2 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2+1) = -A(I1,I2) + B1 * A(I1,I2+1)
                  A(I1,I2) = U
  470          CONTINUE
C
               I1 = J + M1
               IF (I1 .GT. N) GO TO 480
               G = A(I1,1)
               A(I1,1) = B1 * A(I1,1)
  480          IF (.NOT. MATZ) GO TO 500
C
               DO 490 L = 1, N
                  U = B2 * Z(L,J1) + Z(L,J)
                  Z(L,J) = -Z(L,J1) + B1 * Z(L,J)
                  Z(L,J1) = U
  490          CONTINUE
C
  500       CONTINUE
C
  600    CONTINUE
C
         IF (MOD(K,64) .NE. 0) GO TO 700
C     .......... RESCALE TO AVOID UNDERFLOW OR OVERFLOW ..........
         DO 650 J = K, N
            IF (D(J) .GE. DMIN) GO TO 650
            MAXL = MAX(1,MB+1-J)
C
            DO 610 L = MAXL, M1
  610       A(J,L) = DMINRT * A(J,L)
C
            IF (J .EQ. N) GO TO 630
            MAXL = MIN(M1,N-J)
C
            DO 620 L = 1, MAXL
               I1 = J + L
               I2 = MB - L
               A(I1,I2) = DMINRT * A(I1,I2)
  620       CONTINUE
C
  630       IF (.NOT. MATZ) GO TO 645
C
            DO 640 L = 1, N
  640       Z(L,J) = DMINRT * Z(L,J)
C
  645       A(J,MB) = DMIN * A(J,MB)
            D(J) = D(J) / DMIN
  650    CONTINUE
C
  700 CONTINUE
C     .......... FORM SQUARE ROOT OF SCALING MATRIX ..........
  800 DO 810 J = 2, N
  810 E(J) = SQRT(D(J))
C
      IF (.NOT. MATZ) GO TO 840
C
      DO 830 J = 1, N
C
         DO 820 K = 2, N
  820    Z(J,K) = E(K) * Z(J,K)
C
  830 CONTINUE
C
  840 U = 1.0E0
C
      DO 850 J = 2, N
         A(J,M1) = U * E(J) * A(J,M1)
         U = E(J)
         E2(J) = A(J,M1) ** 2
         A(J,MB) = D(J) * A(J,MB)
         D(J) = A(J,MB)
         E(J) = A(J,M1)
  850 CONTINUE
C
      D(1) = A(1,MB)
      E(1) = 0.0E0
      E2(1) = 0.0E0
      GO TO 1001
C
  900 DO 950 J = 1, N
         D(J) = A(J,MB)
         E(J) = 0.0E0
         E2(J) = 0.0E0
  950 CONTINUE
C
 1001 RETURN
      END
