*DECK BQR
      SUBROUTINE BQR (NM, N, MB, A, T, R, IERR, NV, RV)
C***BEGIN PROLOGUE  BQR
C***PURPOSE  Compute some of the eigenvalues of a real symmetric
C            matrix using the QR method with shifts of origin.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A6
C***TYPE      SINGLE PRECISION (BQR-S)
C***KEYWORDS  EIGENVALUES, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BQR,
C     NUM. MATH. 16, 85-92(1970) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 266-272(1971).
C
C     This subroutine finds the eigenvalue of smallest (usually)
C     magnitude of a REAL SYMMETRIC BAND matrix using the
C     QR algorithm with shifts of origin.  Consecutive calls
C     can be made to find further eigenvalues.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, A, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        MB is the (half) band width of the matrix, defined as the
C          number of adjacent diagonals, including the principal
C          diagonal, required to specify the non-zero portion of the
C          lower triangle of the matrix.  MB is an INTEGER variable.
C          MB must be less than or equal to N on first call.
C
C        A contains the lower triangle of the symmetric band input
C          matrix stored as an N by MB array.  Its lowest subdiagonal
C          is stored in the last N+1-MB positions of the first column,
C          its next subdiagonal in the last N+2-MB positions of the
C          second column, further subdiagonals similarly, and finally
C          its principal diagonal in the N positions of the last column.
C          Contents of storages not part of the matrix are arbitrary.
C          On a subsequent call, its output contents from the previous
C          call should be passed.  A is a two-dimensional REAL array,
C          dimensioned A(NM,MB).
C
C        T specifies the shift (of eigenvalues) applied to the diagonal
C          of A in forming the input matrix. What is actually determined
C          is the eigenvalue of A+TI (I is the identity matrix) nearest
C          to T.  On a subsequent call, the output value of T from the
C          previous call should be passed if the next nearest eigenvalue
C          is sought.  T is a REAL variable.
C
C        R should be specified as zero on the first call, and as its
C          output value from the previous call on a subsequent call.
C          It is used to determine when the last row and column of
C          the transformed band matrix can be regarded as negligible.
C          R is a REAL variable.
C
C        NV must be set to the dimension of the array parameter RV
C          as declared in the calling program dimension statement.
C          NV is an INTEGER variable.
C
C     On OUTPUT
C
C        A contains the transformed band matrix.  The matrix A+TI
C          derived from the output parameters is similar to the
C          input A+TI to within rounding errors.  Its last row and
C          column are null (if IERR is zero).
C
C        T contains the computed eigenvalue of A+TI (if IERR is zero),
C          where I is the identity matrix.
C
C        R contains the maximum of its input value and the norm of the
C          last column of the input matrix A.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after a total of 30 iterations.
C
C        RV is a one-dimensional REAL array of dimension NV which is
C          at least (2*MB**2+4*MB-3), used for temporary storage.  The
C          first (3*MB-2) locations correspond to the ALGOL array B,
C          the next (2*MB-1) locations correspond to the ALGOL array H,
C          and the final (2*MB**2-MB) locations correspond to the MB
C          by (2*MB-1) ALGOL array U.
C
C     NOTE. For a subsequent call, N should be replaced by N-1, but
C     MB should not be altered even when it exceeds the current N.
C
C     Calls PYTHAG(A,B) for SQRT(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  BQR
C
      INTEGER I,J,K,L,M,N,II,IK,JK,JM,KJ,KK,KM,LL,MB,MK,MN,MZ
      INTEGER M1,M2,M3,M4,NI,NM,NV,ITS,KJ1,M21,M31,IERR,IMULT
      REAL A(NM,*),RV(*)
      REAL F,G,Q,R,S,T,SCALE
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  BQR
      IERR = 0
      M1 = MIN(MB,N)
      M = M1 - 1
      M2 = M + M
      M21 = M2 + 1
      M3 = M21 + M
      M31 = M3 + 1
      M4 = M31 + M2
      MN = M + N
      MZ = MB - M1
      ITS = 0
C     .......... TEST FOR CONVERGENCE ..........
   40 G = A(N,MB)
      IF (M .EQ. 0) GO TO 360
      F = 0.0E0
C
      DO 50 K = 1, M
         MK = K + MZ
         F = F + ABS(A(N,MK))
   50 CONTINUE
C
      IF (ITS .EQ. 0 .AND. F .GT. R) R = F
      IF (R + F .LE. R) GO TO 360
      IF (ITS .EQ. 30) GO TO 1000
      ITS = ITS + 1
C     .......... FORM SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
      IF (F .GT. 0.25E0 * R .AND. ITS .LT. 5) GO TO 90
      F = A(N,MB-1)
      IF (F .EQ. 0.0E0) GO TO 70
      Q = (A(N-1,MB) - G) / (2.0E0 * F)
      S = PYTHAG(Q,1.0E0)
      G = G - F / (Q + SIGN(S,Q))
   70 T = T + G
C
      DO 80 I = 1, N
   80 A(I,MB) = A(I,MB) - G
C
   90 DO 100 K = M31, M4
  100 RV(K) = 0.0E0
C
      DO 350 II = 1, MN
         I = II - M
         NI = N - II
         IF (NI .LT. 0) GO TO 230
C     .......... FORM COLUMN OF SHIFTED MATRIX A-G*I ..........
         L = MAX(1,2-I)
C
         DO 110 K = 1, M3
  110    RV(K) = 0.0E0
C
         DO 120 K = L, M1
            KM = K + M
            MK = K + MZ
            RV(KM) = A(II,MK)
  120    CONTINUE
C
         LL = MIN(M,NI)
         IF (LL .EQ. 0) GO TO 135
C
         DO 130 K = 1, LL
            KM = K + M21
            IK = II + K
            MK = MB - K
            RV(KM) = A(IK,MK)
  130    CONTINUE
C     .......... PRE-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
  135    LL = M2
         IMULT = 0
C     .......... MULTIPLICATION PROCEDURE ..........
  140    KJ = M4 - M1
C
         DO 170 J = 1, LL
            KJ = KJ + M1
            JM = J + M3
            IF (RV(JM) .EQ. 0.0E0) GO TO 170
            F = 0.0E0
C
            DO 150 K = 1, M1
               KJ = KJ + 1
               JK = J + K - 1
               F = F + RV(KJ) * RV(JK)
  150       CONTINUE
C
            F = F / RV(JM)
            KJ = KJ - M1
C
            DO 160 K = 1, M1
               KJ = KJ + 1
               JK = J + K - 1
               RV(JK) = RV(JK) - RV(KJ) * F
  160       CONTINUE
C
            KJ = KJ - M1
  170    CONTINUE
C
         IF (IMULT .NE. 0) GO TO 280
C     .......... HOUSEHOLDER REFLECTION ..........
         F = RV(M21)
         S = 0.0E0
         RV(M4) = 0.0E0
         SCALE = 0.0E0
C
         DO 180 K = M21, M3
  180    SCALE = SCALE + ABS(RV(K))
C
         IF (SCALE .EQ. 0.0E0) GO TO 210
C
         DO 190 K = M21, M3
  190    S = S + (RV(K)/SCALE)**2
C
         S = SCALE * SCALE * S
         G = -SIGN(SQRT(S),F)
         RV(M21) = G
         RV(M4) = S - F * G
         KJ = M4 + M2 * M1 + 1
         RV(KJ) = F - G
C
         DO 200 K = 2, M1
            KJ = KJ + 1
            KM = K + M2
            RV(KJ) = RV(KM)
  200    CONTINUE
C     .......... SAVE COLUMN OF TRIANGULAR FACTOR R ..........
  210    DO 220 K = L, M1
            KM = K + M
            MK = K + MZ
            A(II,MK) = RV(KM)
  220    CONTINUE
C
  230    L = MAX(1,M1+1-I)
         IF (I .LE. 0) GO TO 300
C     .......... PERFORM ADDITIONAL STEPS ..........
         DO 240 K = 1, M21
  240    RV(K) = 0.0E0
C
         LL = MIN(M1,NI+M1)
C     .......... GET ROW OF TRIANGULAR FACTOR R ..........
         DO 250 KK = 1, LL
            K = KK - 1
            KM = K + M1
            IK = I + K
            MK = MB - K
            RV(KM) = A(IK,MK)
  250    CONTINUE
C     .......... POST-MULTIPLY WITH HOUSEHOLDER REFLECTIONS ..........
         LL = M1
         IMULT = 1
         GO TO 140
C     .......... STORE COLUMN OF NEW A MATRIX ..........
  280    DO 290 K = L, M1
            MK = K + MZ
            A(I,MK) = RV(K)
  290    CONTINUE
C     .......... UPDATE HOUSEHOLDER REFLECTIONS ..........
  300    IF (L .GT. 1) L = L - 1
         KJ1 = M4 + L * M1
C
         DO 320 J = L, M2
            JM = J + M3
            RV(JM) = RV(JM+1)
C
            DO 320 K = 1, M1
               KJ1 = KJ1 + 1
               KJ = KJ1 - M1
               RV(KJ) = RV(KJ1)
  320    CONTINUE
C
  350 CONTINUE
C
      GO TO 40
C     .......... CONVERGENCE ..........
  360 T = T + G
C
      DO 380 I = 1, N
  380 A(I,MB) = A(I,MB) - G
C
      DO 400 K = 1, M1
         MK = K + MZ
         A(N,MK) = 0.0E0
  400 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = N
 1001 RETURN
      END
