*DECK MINFIT
      SUBROUTINE MINFIT (NM, M, N, A, W, IP, B, IERR, RV1)
C***BEGIN PROLOGUE  MINFIT
C***PURPOSE  Compute the singular value decomposition of a rectangular
C            matrix and solve the related linear least squares problem.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D9
C***TYPE      SINGLE PRECISION (MINFIT-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure MINFIT,
C     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     This subroutine determines, towards the solution of the linear
C                                                        T
C     system AX=B, the singular value decomposition A=USV  of a real
C                                         T
C     M by N rectangular matrix, forming U B rather than U.  Householder
C     bidiagonalization and a variant of the QR algorithm are used.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and B, as declared in the calling
C          program dimension statement.  Note that NM must be at least
C          as large as the maximum of M and N.  NM is an INTEGER
C          variable.
C
C        M is the number of rows of A and B.  M is an INTEGER variable.
C
C        N is the number of columns of A and the order of V.  N is an
C          INTEGER variable.
C
C        A contains the rectangular coefficient matrix of the system.
C          A is a two-dimensional REAL array, dimensioned A(NM,N).
C
C        IP is the number of columns of B.  IP can be zero.
C
C        B contains the constant column matrix of the system if IP is
C          not zero.  Otherwise, B is not referenced.  B is a two-
C          dimensional REAL array, dimensioned B(NM,IP).
C
C     On OUTPUT
C
C        A has been overwritten by the matrix V (orthogonal) of the
C          decomposition in its first N rows and columns.  If an
C          error exit is made, the columns of V corresponding to
C          indices of correct singular values should be correct.
C
C        W contains the N (non-negative) singular values of A (the
C          diagonal elements of S).  They are unordered.  If an
C          error exit is made, the singular values should be correct
C          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
C          REAL array, dimensioned W(N).
C
C                                   T
C        B has been overwritten by U B.  If an error exit is made,
C                       T
C          the rows of U B corresponding to indices of correct singular
C          values should be correct.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          K          if the K-th singular value has not been
C                     determined after 30 iterations.
C                     The singular values should be correct for
C                     indices IERR+1, IERR+2, ..., N.
C
C        RV1 is a one-dimensional REAL array used for temporary storage,
C          dimensioned RV1(N).
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
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
C***END PROLOGUE  MINFIT
C
      INTEGER I,J,K,L,M,N,II,IP,I1,KK,K1,LL,L1,M1,NM,ITS,IERR
      REAL A(NM,*),W(*),B(NM,IP),RV1(*)
      REAL C,F,G,H,S,X,Y,Z,SCALE,S1
      REAL PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  MINFIT
      IERR = 0
C     .......... HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM ..........
      G = 0.0E0
      SCALE = 0.0E0
      S1 = 0.0E0
C
      DO 300 I = 1, N
         L = I + 1
         RV1(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M) GO TO 210
C
         DO 120 K = I, M
  120    SCALE = SCALE + ABS(A(K,I))
C
         IF (SCALE .EQ. 0.0E0) GO TO 210
C
         DO 130 K = I, M
            A(K,I) = A(K,I) / SCALE
            S = S + A(K,I)**2
  130    CONTINUE
C
         F = A(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         A(I,I) = F - G
         IF (I .EQ. N) GO TO 160
C
         DO 150 J = L, N
            S = 0.0E0
C
            DO 140 K = I, M
  140       S = S + A(K,I) * A(K,J)
C
            F = S / H
C
            DO 150 K = I, M
               A(K,J) = A(K,J) + F * A(K,I)
  150    CONTINUE
C
  160    IF (IP .EQ. 0) GO TO 190
C
         DO 180 J = 1, IP
            S = 0.0E0
C
            DO 170 K = I, M
  170       S = S + A(K,I) * B(K,J)
C
            F = S / H
C
            DO 180 K = I, M
               B(K,J) = B(K,J) + F * A(K,I)
  180    CONTINUE
C
  190    DO 200 K = I, M
  200    A(K,I) = SCALE * A(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N
  220    SCALE = SCALE + ABS(A(I,K))
C
         IF (SCALE .EQ. 0.0E0) GO TO 290
C
         DO 230 K = L, N
            A(I,K) = A(I,K) / SCALE
            S = S + A(I,K)**2
  230    CONTINUE
C
         F = A(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         A(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = A(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M
            S = 0.0E0
C
            DO 250 K = L, N
  250       S = S + A(J,K) * A(I,K)
C
            DO 260 K = L, N
               A(J,K) = A(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N
  280    A(I,K) = SCALE * A(I,K)
C
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
C     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
C                FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0E0) GO TO 360
C
         DO 320 J = L, N
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    A(J,I) = (A(I,J) / A(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0E0
C
            DO 340 K = L, N
  340       S = S + A(I,K) * A(K,J)
C
            DO 350 K = L, N
               A(K,J) = A(K,J) + S * A(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            A(I,J) = 0.0E0
            A(J,I) = 0.0E0
  380    CONTINUE
C
  390    A(I,I) = 1.0E0
         G = RV1(I)
         L = I
  400 CONTINUE
C
      IF (M .GE. N .OR. IP .EQ. 0) GO TO 510
      M1 = M + 1
C
      DO 500 I = M1, N
C
         DO 500 J = 1, IP
            B(I,J) = 0.0E0
  500 CONTINUE
C     .......... DIAGONALIZATION OF THE BIDIAGONAL FORM ..........
  510 CONTINUE
C     .......... FOR K=N STEP -1 UNTIL 1 DO -- ..........
      DO 700 KK = 1, N
         K1 = N - KK
         K = K1 + 1
         ITS = 0
C     .......... TEST FOR SPLITTING.
C                FOR L=K STEP -1 UNTIL 1 DO -- ..........
  520    DO 530 LL = 1, K
            L1 = K - LL
            L = L1 + 1
            IF (S1 + ABS(RV1(L)) .EQ. S1) GO TO 565
C     .......... RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
            IF (S1 + ABS(W(L1)) .EQ. S1) GO TO 540
  530    CONTINUE
C     .......... CANCELLATION OF RV1(L) IF L GREATER THAN 1 ..........
  540    C = 0.0E0
         S = 1.0E0
C
         DO 560 I = L, K
            F = S * RV1(I)
            RV1(I) = C * RV1(I)
            IF (S1 + ABS(F) .EQ. S1) GO TO 565
            G = W(I)
            H = PYTHAG(F,G)
            W(I) = H
            C = G / H
            S = -F / H
            IF (IP .EQ. 0) GO TO 560
C
            DO 550 J = 1, IP
               Y = B(L1,J)
               Z = B(I,J)
               B(L1,J) = Y * C + Z * S
               B(I,J) = -Y * S + Z * C
  550       CONTINUE
C
  560    CONTINUE
C     .......... TEST FOR CONVERGENCE ..........
  565    Z = W(K)
         IF (L .EQ. K) GO TO 650
C     .......... SHIFT FROM BOTTOM 2 BY 2 MINOR ..........
         IF (ITS .EQ. 30) GO TO 1000
         ITS = ITS + 1
         X = W(L)
         Y = W(K1)
         G = RV1(K1)
         H = RV1(K)
         F = 0.5E0 * (((G + Z) / H) * ((G - Z) / Y) + Y / H - H / Y)
         G = PYTHAG(F,1.0E0)
         F = X - (Z / X) * Z + (H / X) * (Y / (F + SIGN(G,F)) - H)
C     .......... NEXT QR TRANSFORMATION ..........
         C = 1.0E0
         S = 1.0E0
C
         DO 600 I1 = L, K1
            I = I1 + 1
            G = RV1(I)
            Y = W(I)
            H = S * G
            G = C * G
            Z = PYTHAG(F,H)
            RV1(I1) = Z
            C = F / Z
            S = H / Z
            F = X * C + G * S
            G = -X * S + G * C
            H = Y * S
            Y = Y * C
C
            DO 570 J = 1, N
               X = A(J,I1)
               Z = A(J,I)
               A(J,I1) = X * C + Z * S
               A(J,I) = -X * S + Z * C
  570       CONTINUE
C
            Z = PYTHAG(F,H)
            W(I1) = Z
C     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0E0) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (IP .EQ. 0) GO TO 600
C
            DO 590 J = 1, IP
               Y = B(I1,J)
               Z = B(I,J)
               B(I1,J) = Y * C + Z * S
               B(I,J) = -Y * S + Z * C
  590       CONTINUE
C
  600    CONTINUE
C
         RV1(L) = 0.0E0
         RV1(K) = F
         W(K) = X
         GO TO 520
C     .......... CONVERGENCE ..........
  650    IF (Z .GE. 0.0E0) GO TO 700
C     .......... W(K) IS MADE NON-NEGATIVE ..........
         W(K) = -Z
C
         DO 690 J = 1, N
  690    A(J,K) = -A(J,K)
C
  700 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO A
C                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END
