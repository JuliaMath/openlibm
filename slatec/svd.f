*DECK SVD
      SUBROUTINE SVD (NM, M, N, A, W, MATU, U, MATV, V, IERR, RV1)
C***BEGIN PROLOGUE  SVD
C***SUBSIDIARY
C***PURPOSE  Perform the singular value decomposition of a rectangular
C            matrix.
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (SVD-S)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure SVD,
C     NUM. MATH. 14, 403-420(1970) by Golub and Reinsch.
C     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
C
C     This subroutine determines the singular value decomposition
C          T
C     A=USV  of a REAL M by N rectangular matrix.  Householder
C     bidiagonalization and a variant of the QR algorithm are used.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A, U and V, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C          Note that NM must be at least as large as the maximum
C          of M and N.
C
C        M is the number of rows of A and U.
C
C        N is the number of columns of A and U and the order of V.
C
C        A contains the rectangular input matrix to be decomposed.  A is
C          a two-dimensional REAL array, dimensioned A(NM,N).
C
C        MATU should be set to .TRUE. if the U matrix in the
C          decomposition is desired, and to .FALSE. otherwise.
C          MATU is a LOGICAL variable.
C
C        MATV should be set to .TRUE. if the V matrix in the
C          decomposition is desired, and to .FALSE. otherwise.
C          MATV is a LOGICAL variable.
C
C     On Output
C
C        A is unaltered (unless overwritten by U or V).
C
C        W contains the N (non-negative) singular values of A (the
C          diagonal elements of S).  They are unordered.  If an
C          error exit is made, the singular values should be correct
C          for indices IERR+1, IERR+2, ..., N.  W is a one-dimensional
C          REAL array, dimensioned W(N).
C
C        U contains the matrix U (orthogonal column vectors) of the
C          decomposition if MATU has been set to .TRUE.  Otherwise,
C          U is used as a temporary array.  U may coincide with A.
C          If an error exit is made, the columns of U corresponding
C          to indices of correct singular values should be correct.
C          U is a two-dimensional REAL array, dimensioned U(NM,N).
C
C        V contains the matrix V (orthogonal) of the decomposition if
C          MATV has been set to .TRUE.  Otherwise, V is not referenced.
C          V may also coincide with A if U does not.  If an error
C          exit is made, the columns of V corresponding to indices of
C          correct singular values should be correct.  V is a two-
C          dimensional REAL array, dimensioned V(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          K          if the K-th singular value has not been
C                     determined after 30 iterations.
C
C        RV1 is a one-dimensional REAL array used for temporary storage,
C          dimensioned RV1(N).
C
C     CALLS PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***SEE ALSO  EISDOC
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   811101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  SVD
C
      INTEGER I,J,K,L,M,N,II,I1,KK,K1,LL,L1,MN,NM,ITS,IERR
      REAL A(NM,*),W(*),U(NM,*),V(NM,*),RV1(*)
      REAL C,F,G,H,S,X,Y,Z,SCALE,S1
      REAL PYTHAG
      LOGICAL MATU,MATV
C
C***FIRST EXECUTABLE STATEMENT  SVD
      IERR = 0
C
      DO 100 I = 1, M
C
         DO 100 J = 1, N
            U(I,J) = A(I,J)
  100 CONTINUE
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
  120    SCALE = SCALE + ABS(U(K,I))
C
         IF (SCALE .EQ. 0.0E0) GO TO 210
C
         DO 130 K = I, M
            U(K,I) = U(K,I) / SCALE
            S = S + U(K,I)**2
  130    CONTINUE
C
         F = U(I,I)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,I) = F - G
         IF (I .EQ. N) GO TO 190
C
         DO 150 J = L, N
            S = 0.0E0
C
            DO 140 K = I, M
  140       S = S + U(K,I) * U(K,J)
C
            F = S / H
C
            DO 150 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  150    CONTINUE
C
  190    DO 200 K = I, M
  200    U(K,I) = SCALE * U(K,I)
C
  210    W(I) = SCALE * G
         G = 0.0E0
         S = 0.0E0
         SCALE = 0.0E0
         IF (I .GT. M .OR. I .EQ. N) GO TO 290
C
         DO 220 K = L, N
  220    SCALE = SCALE + ABS(U(I,K))
C
         IF (SCALE .EQ. 0.0E0) GO TO 290
C
         DO 230 K = L, N
            U(I,K) = U(I,K) / SCALE
            S = S + U(I,K)**2
  230    CONTINUE
C
         F = U(I,L)
         G = -SIGN(SQRT(S),F)
         H = F * G - S
         U(I,L) = F - G
C
         DO 240 K = L, N
  240    RV1(K) = U(I,K) / H
C
         IF (I .EQ. M) GO TO 270
C
         DO 260 J = L, M
            S = 0.0E0
C
            DO 250 K = L, N
  250       S = S + U(J,K) * U(I,K)
C
            DO 260 K = L, N
               U(J,K) = U(J,K) + S * RV1(K)
  260    CONTINUE
C
  270    DO 280 K = L, N
  280    U(I,K) = SCALE * U(I,K)
C
  290    S1 = MAX(S1,ABS(W(I))+ABS(RV1(I)))
  300 CONTINUE
C     .......... ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS ..........
      IF (.NOT. MATV) GO TO 410
C     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO 400 II = 1, N
         I = N + 1 - II
         IF (I .EQ. N) GO TO 390
         IF (G .EQ. 0.0E0) GO TO 360
C
         DO 320 J = L, N
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
  320    V(J,I) = (U(I,J) / U(I,L)) / G
C
         DO 350 J = L, N
            S = 0.0E0
C
            DO 340 K = L, N
  340       S = S + U(I,K) * V(K,J)
C
            DO 350 K = L, N
               V(K,J) = V(K,J) + S * V(K,I)
  350    CONTINUE
C
  360    DO 380 J = L, N
            V(I,J) = 0.0E0
            V(J,I) = 0.0E0
  380    CONTINUE
C
  390    V(I,I) = 1.0E0
         G = RV1(I)
         L = I
  400 CONTINUE
C     .......... ACCUMULATION OF LEFT-HAND TRANSFORMATIONS ..........
  410 IF (.NOT. MATU) GO TO 510
C     ..........FOR I=MIN(M,N) STEP -1 UNTIL 1 DO -- ..........
      MN = N
      IF (M .LT. N) MN = M
C
      DO 500 II = 1, MN
         I = MN + 1 - II
         L = I + 1
         G = W(I)
         IF (I .EQ. N) GO TO 430
C
         DO 420 J = L, N
  420    U(I,J) = 0.0E0
C
  430    IF (G .EQ. 0.0E0) GO TO 475
         IF (I .EQ. MN) GO TO 460
C
         DO 450 J = L, N
            S = 0.0E0
C
            DO 440 K = L, M
  440       S = S + U(K,I) * U(K,J)
C     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            F = (S / U(I,I)) / G
C
            DO 450 K = I, M
               U(K,J) = U(K,J) + F * U(K,I)
  450    CONTINUE
C
  460    DO 470 J = I, M
  470    U(J,I) = U(J,I) / G
C
         GO TO 490
C
  475    DO 480 J = I, M
  480    U(J,I) = 0.0E0
C
  490    U(I,I) = U(I,I) + 1.0E0
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
            IF (.NOT. MATU) GO TO 560
C
            DO 550 J = 1, M
               Y = U(J,L1)
               Z = U(J,I)
               U(J,L1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
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
            IF (.NOT. MATV) GO TO 575
C
            DO 570 J = 1, N
               X = V(J,I1)
               Z = V(J,I)
               V(J,I1) = X * C + Z * S
               V(J,I) = -X * S + Z * C
  570       CONTINUE
C
  575       Z = PYTHAG(F,H)
            W(I1) = Z
C     .......... ROTATION CAN BE ARBITRARY IF Z IS ZERO ..........
            IF (Z .EQ. 0.0E0) GO TO 580
            C = F / Z
            S = H / Z
  580       F = C * G + S * Y
            X = -S * G + C * Y
            IF (.NOT. MATU) GO TO 600
C
            DO 590 J = 1, M
               Y = U(J,I1)
               Z = U(J,I)
               U(J,I1) = Y * C + Z * S
               U(J,I) = -Y * S + Z * C
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
         IF (.NOT. MATV) GO TO 700
C
         DO 690 J = 1, N
  690    V(J,K) = -V(J,K)
C
  700 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO A
C                SINGULAR VALUE AFTER 30 ITERATIONS ..........
 1000 IERR = K
 1001 RETURN
      END
