*DECK LSI
      SUBROUTINE LSI (W, MDW, MA, MG, N, PRGOPT, X, RNORM, MODE, WS, IP)
C***BEGIN PROLOGUE  LSI
C***SUBSIDIARY
C***PURPOSE  Subsidiary to LSEI
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LSI-S, DLSI-D)
C***AUTHOR  Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     This is a companion subprogram to LSEI.  The documentation for
C     LSEI has complete usage instructions.
C
C     Solve..
C              AX = B,  A  MA by N  (least squares equations)
C     subject to..
C
C              GX.GE.H, G  MG by N  (inequality constraints)
C
C     Input..
C
C      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
C                       (G H)
C
C     MDW,MA,MG,N
C              contain (resp) var. dimension of W(*,*),
C              and matrix dimensions.
C
C     PRGOPT(*),
C              Program option vector.
C
C     OUTPUT..
C
C      X(*),RNORM
C
C              Solution vector(unless MODE=2), length of AX-B.
C
C      MODE
C              =0   Inequality constraints are compatible.
C              =2   Inequality constraints contradictory.
C
C      WS(*),
C              Working storage of dimension K+N+(MG+2)*(N+7),
C              where K=MAX(MA+MG,N).
C      IP(MG+2*N+1)
C              Integer working storage
C
C***ROUTINES CALLED  H12, HFTI, LPDP, R1MACH, SASUM, SAXPY, SCOPY, SDOT,
C                    SSCAL, SSWAP
C***REVISION HISTORY  (YYMMDD)
C   790701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890618  Completely restructured and extensively revised (WRB & RWC)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   920422  Changed CALL to HFTI to include variable MA.  (WRB)
C***END PROLOGUE  LSI
      INTEGER IP(*), MA, MDW, MG, MODE, N
      REAL             PRGOPT(*), RNORM, W(MDW,*), WS(*), X(*)
C
      EXTERNAL H12, HFTI, LPDP, R1MACH, SASUM, SAXPY, SCOPY, SDOT,
     *   SSCAL, SSWAP
      REAL             R1MACH, SASUM, SDOT
C
      REAL             ANORM, FAC, GAM, RB, SRELPR, TAU, TOL, XNORM
      INTEGER I, J, K, KEY, KRANK, KRM1, KRP1, L, LAST, LINK, M, MAP1,
     *   MDLPDP, MINMAN, N1, N2, N3, NEXT, NP1
      LOGICAL COV, FIRST, SCLCOV
C
      SAVE SRELPR, FIRST
      DATA FIRST /.TRUE./
C
C***FIRST EXECUTABLE STATEMENT  LSI
C
C     Set the nominal tolerance used in the code.
C
      IF (FIRST) SRELPR = R1MACH(4)
      FIRST = .FALSE.
      TOL = SQRT(SRELPR)
C
      MODE = 0
      RNORM = 0.E0
      M = MA + MG
      NP1 = N + 1
      KRANK = 0
      IF (N.LE.0 .OR. M.LE.0) GO TO 370
C
C     To process option vector.
C
      COV = .FALSE.
      SCLCOV = .TRUE.
      LAST = 1
      LINK = PRGOPT(1)
C
  100 IF (LINK.GT.1) THEN
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) COV = PRGOPT(LAST+2) .NE. 0.E0
         IF (KEY.EQ.10) SCLCOV = PRGOPT(LAST+2) .EQ. 0.E0
         IF (KEY.EQ.5) TOL = MAX(SRELPR,PRGOPT(LAST+2))
         NEXT = PRGOPT(LINK)
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
C
C     Compute matrix norm of least squares equations.
C
      ANORM = 0.E0
      DO 110 J = 1,N
         ANORM = MAX(ANORM,SASUM(MA,W(1,J),1))
  110 CONTINUE
C
C     Set tolerance for HFTI( ) rank test.
C
      TAU = TOL*ANORM
C
C     Compute Householder orthogonal decomposition of matrix.
C
      CALL SCOPY (N, 0.E0, 0, WS, 1)
      CALL SCOPY (MA, W(1, NP1), 1, WS, 1)
      K = MAX(M,N)
      MINMAN = MIN(MA,N)
      N1 = K + 1
      N2 = N1 + N
      CALL HFTI (W, MDW, MA, N, WS, MA, 1, TAU, KRANK, RNORM, WS(N2),
     +           WS(N1), IP)
      FAC = 1.E0
      GAM = MA - KRANK
      IF (KRANK.LT.MA .AND. SCLCOV) FAC = RNORM**2/GAM
C
C     Reduce to LPDP and solve.
C
      MAP1 = MA + 1
C
C     Compute inequality rt-hand side for LPDP.
C
      IF (MA.LT.M) THEN
         IF (MINMAN.GT.0) THEN
            DO 120 I = MAP1,M
               W(I,NP1) = W(I,NP1) - SDOT(N,W(I,1),MDW,WS,1)
  120       CONTINUE
C
C           Apply permutations to col. of inequality constraint matrix.
C
            DO 130 I = 1,MINMAN
               CALL SSWAP (MG, W(MAP1,I), 1, W(MAP1,IP(I)), 1)
  130       CONTINUE
C
C           Apply Householder transformations to constraint matrix.
C
            IF (KRANK.GT.0 .AND. KRANK.LT.N) THEN
               DO 140 I = KRANK,1,-1
                  CALL H12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),
     +                      W(MAP1,1), MDW, 1, MG)
  140          CONTINUE
            ENDIF
C
C           Compute permuted inequality constraint matrix times r-inv.
C
            DO 160 I = MAP1,M
               DO 150 J = 1,KRANK
                  W(I,J) = (W(I,J)-SDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  150          CONTINUE
  160       CONTINUE
         ENDIF
C
C        Solve the reduced problem with LPDP algorithm,
C        the least projected distance problem.
C
         CALL LPDP(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X,
     +             XNORM, MDLPDP, WS(N2), IP(N+1))
C
C        Compute solution in original coordinates.
C
         IF (MDLPDP.EQ.1) THEN
            DO 170 I = KRANK,1,-1
               X(I) = (X(I)-SDOT(KRANK-I,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170       CONTINUE
C
C           Apply Householder transformation to solution vector.
C
            IF (KRANK.LT.N) THEN
               DO 180 I = 1,KRANK
                  CALL H12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),
     +                      X, 1, 1, 1)
  180          CONTINUE
            ENDIF
C
C           Repermute variables to their input order.
C
            IF (MINMAN.GT.0) THEN
               DO 190 I = MINMAN,1,-1
                  CALL SSWAP (1, X(I), 1, X(IP(I)), 1)
  190          CONTINUE
C
C              Variables are now in original coordinates.
C              Add solution of unconstrained problem.
C
               DO 200 I = 1,N
                  X(I) = X(I) + WS(I)
  200          CONTINUE
C
C              Compute the residual vector norm.
C
               RNORM = SQRT(RNORM**2+XNORM**2)
            ENDIF
         ELSE
            MODE = 2
         ENDIF
      ELSE
         CALL SCOPY (N, WS, 1, X, 1)
      ENDIF
C
C     Compute covariance matrix based on the orthogonal decomposition
C     from HFTI( ).
C
      IF (.NOT.COV .OR. KRANK.LE.0) GO TO 370
      KRM1 = KRANK - 1
      KRP1 = KRANK + 1
C
C     Copy diagonal terms to working array.
C
      CALL SCOPY (KRANK, W, MDW+1, WS(N2), 1)
C
C     Reciprocate diagonal terms.
C
      DO 210 J = 1,KRANK
         W(J,J) = 1.E0/W(J,J)
  210 CONTINUE
C
C     Invert the upper triangular QR factor on itself.
C
      IF (KRANK.GT.1) THEN
         DO 230 I = 1,KRM1
            DO 220 J = I+1,KRANK
               W(I,J) = -SDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  220       CONTINUE
  230    CONTINUE
      ENDIF
C
C     Compute the inverted factor times its transpose.
C
      DO 250 I = 1,KRANK
         DO 240 J = I,KRANK
            W(I,J) = SDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  240    CONTINUE
  250 CONTINUE
C
C     Zero out lower trapezoidal part.
C     Copy upper triangular to lower triangular part.
C
      IF (KRANK.LT.N) THEN
         DO 260 J = 1,KRANK
            CALL SCOPY (J, W(1,J), 1, W(J,1), MDW)
  260    CONTINUE
C
         DO 270 I = KRP1,N
            CALL SCOPY (I, 0.E0, 0, W(I,1), MDW)
  270    CONTINUE
C
C        Apply right side transformations to lower triangle.
C
         N3 = N2 + KRP1
         DO 330 I = 1,KRANK
            L = N1 + I
            K = N2 + I
            RB = WS(L-1)*WS(K-1)
C
C           If RB.GE.0.E0, transformation can be regarded as zero.
C
            IF (RB.LT.0.E0) THEN
               RB = 1.E0/RB
C
C              Store unscaled rank one Householder update in work array.
C
               CALL SCOPY (N, 0.E0, 0, WS(N3), 1)
               L = N1 + I
               K = N3 + I
               WS(K-1) = WS(L-1)
C
               DO 280 J = KRP1,N
                  WS(N3+J-1) = W(I,J)
  280          CONTINUE
C
               DO 290 J = 1,N
                  WS(J) = RB*(SDOT(J-I,W(J,I),MDW,WS(N3+I-1),1)+
     +                    SDOT(N-J+1,W(J,J),1,WS(N3+J-1),1))
  290          CONTINUE
C
               L = N3 + I
               GAM = 0.5E0*RB*SDOT(N-I+1,WS(L-1),1,WS(I),1)
               CALL SAXPY (N-I+1, GAM, WS(L-1), 1, WS(I), 1)
               DO 320 J = I,N
                  DO 300 L = 1,I-1
                     W(J,L) = W(J,L) + WS(N3+J-1)*WS(L)
  300             CONTINUE
C
                  DO 310 L = I,J
                     W(J,L) = W(J,L) + WS(J)*WS(N3+L-1)+WS(L)*WS(N3+J-1)
  310             CONTINUE
  320          CONTINUE
            ENDIF
  330    CONTINUE
C
C        Copy lower triangle to upper triangle to symmetrize the
C        covariance matrix.
C
         DO 340 I = 1,N
            CALL SCOPY (I, W(I,1), MDW, W(1,I), 1)
  340    CONTINUE
      ENDIF
C
C     Repermute rows and columns.
C
      DO 350 I = MINMAN,1,-1
         K = IP(I)
         IF (I.NE.K) THEN
            CALL SSWAP (1, W(I,I), 1, W(K,K), 1)
            CALL SSWAP (I-1, W(1,I), 1, W(1,K), 1)
            CALL SSWAP (K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
            CALL SSWAP (N-K, W(I, K+1), MDW, W(K, K+1), MDW)
         ENDIF
  350 CONTINUE
C
C     Put in normalized residual sum of squares scale factor
C     and symmetrize the resulting covariance matrix.
C
      DO 360 J = 1,N
         CALL SSCAL (J, FAC, W(1,J), 1)
         CALL SCOPY (J, W(1,J), 1, W(J,1), MDW)
  360 CONTINUE
C
  370 IP(1) = KRANK
      IP(2) = N + MAX(M,N) + (MG+2)*(N+7)
      RETURN
      END
