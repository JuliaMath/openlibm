*DECK CINVIT
      SUBROUTINE CINVIT (NM, N, AR, AI, WR, WI, SELECT, MM, M, ZR, ZI,
     +   IERR, RM1, RM2, RV1, RV2)
C***BEGIN PROLOGUE  CINVIT
C***PURPOSE  Compute the eigenvectors of a complex upper Hessenberg
C            associated with specified eigenvalues using inverse
C            iteration.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      COMPLEX (INVIT-S, CINVIT-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure CXINVIT
C     by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP. VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     This subroutine finds those eigenvectors of A COMPLEX UPPER
C     Hessenberg matrix corresponding to specified eigenvalues,
C     using inverse iteration.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, AR, AI, ZR and ZI, as declared in the
C          calling program dimension statement.  NM is an INTEGER
C          variable.
C
C        N is the order of the matrix A=(AR,AI).  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        AR and AI contain the real and imaginary parts, respectively,
C          of the complex upper Hessenberg matrix.  AR and AI are
C          two-dimensional REAL arrays, dimensioned AR(NM,N)
C          and AI(NM,N).
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues of the matrix.  The eigenvalues must be
C          stored in a manner identical to that of subroutine  COMLR,
C          which recognizes possible splitting of the matrix.  WR and
C          WI are one-dimensional REAL arrays, dimensioned WR(N) and
C          WI(N).
C
C        SELECT specifies the eigenvectors to be found.  The
C          eigenvector corresponding to the J-th eigenvalue is
C          specified by setting SELECT(J) to .TRUE.  SELECT is a
C          one-dimensional LOGICAL array, dimensioned SELECT(N).
C
C        MM should be set to an upper bound for the number of
C          eigenvectors to be found.  MM is an INTEGER variable.
C
C     On OUTPUT
C
C        AR, AI, WI, and SELECT are unaltered.
C
C        WR may have been altered since close eigenvalues are perturbed
C          slightly in searching for independent eigenvectors.
C
C        M is the number of eigenvectors actually found.  M is an
C          INTEGER variable.
C
C        ZR and ZI contain the real and imaginary parts, respectively,
C          of the eigenvectors corresponding to the flagged eigenvalues.
C          The eigenvectors are normalized so that the component of
C          largest magnitude is 1.  Any vector which fails the
C          acceptance test is set to zero.  ZR and ZI are
C          two-dimensional REAL arrays, dimensioned ZR(NM,MM) and
C          ZI(NM,MM).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          -(2*N+1)   if more than MM eigenvectors have been requested
C                     (the MM eigenvectors calculated to this point are
C                     in ZR and ZI),
C          -K         if the iteration corresponding to the K-th
C                     value fails (if this occurs more than once, K
C                     is the index of the last occurrence); the
C                     corresponding columns of ZR and ZI are set to
C                     zero vectors,
C          -(N+K)     if both error situations occur.
C
C        RV1 and RV2 are one-dimensional REAL arrays used for
C          temporary storage, dimensioned RV1(N) and RV2(N).
C          They hold the approximate eigenvectors during the inverse
C          iteration process.
C
C        RM1 and RM2 are two-dimensional REAL arrays used for
C          temporary storage, dimensioned RM1(N,N) and RM2(N,N).
C          These arrays hold the triangularized form of the upper
C          Hessenberg matrix used in the inverse iteration process.
C
C     The ALGOL procedure GUESSVEC appears in CINVIT in-line.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C     Calls CDIV for complex division.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  CDIV, PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CINVIT
C
      INTEGER I,J,K,M,N,S,II,MM,MP,NM,UK,IP1,ITS,KM1,IERR
      REAL AR(NM,*),AI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
      REAL RM1(N,*),RM2(N,*),RV1(*),RV2(*)
      REAL X,Y,EPS3,NORM,NORMV,GROWTO,ILAMBD,RLAMBD,UKROOT
      REAL PYTHAG
      LOGICAL SELECT(N)
C
C***FIRST EXECUTABLE STATEMENT  CINVIT
      IERR = 0
      UK = 0
      S = 1
C
      DO 980 K = 1, N
         IF (.NOT. SELECT(K)) GO TO 980
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (AR(UK+1,UK) .EQ. 0.0E0 .AND. AI(UK+1,UK) .EQ. 0.0E0)
     1         GO TO 140
  120    CONTINUE
C     .......... COMPUTE INFINITY NORM OF LEADING UK BY UK
C                (HESSENBERG) MATRIX ..........
  140    NORM = 0.0E0
         MP = 1
C
         DO 180 I = 1, UK
            X = 0.0E0
C
            DO 160 J = MP, UK
  160       X = X + PYTHAG(AR(I,J),AI(I,J))
C
            IF (X .GT. NORM) NORM = X
            MP = I
  180    CONTINUE
C     .......... EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION
C                AND CLOSE ROOTS ARE MODIFIED BY EPS3 ..........
         IF (NORM .EQ. 0.0E0) NORM = 1.0E0
         EPS3 = NORM
  190    EPS3 = 0.5E0*EPS3
         IF (NORM + EPS3 .GT. NORM) GO TO 190
         EPS3 = 2.0E0*EPS3
C     .......... GROWTO IS THE CRITERION FOR GROWTH ..........
         UKROOT = SQRT(REAL(UK))
         GROWTO = 0.1E0 / UKROOT
  200    RLAMBD = WR(K)
         ILAMBD = WI(K)
         IF (K .EQ. 1) GO TO 280
         KM1 = K - 1
         GO TO 240
C     .......... PERTURB EIGENVALUE IF IT IS CLOSE
C                TO ANY PREVIOUS EIGENVALUE ..........
  220    RLAMBD = RLAMBD + EPS3
C     .......... FOR I=K-1 STEP -1 UNTIL 1 DO -- ..........
  240    DO 260 II = 1, KM1
            I = K - II
            IF (SELECT(I) .AND. ABS(WR(I)-RLAMBD) .LT. EPS3 .AND.
     1         ABS(WI(I)-ILAMBD) .LT. EPS3) GO TO 220
  260    CONTINUE
C
         WR(K) = RLAMBD
C     .......... FORM UPPER HESSENBERG (AR,AI)-(RLAMBD,ILAMBD)*I
C                AND INITIAL COMPLEX VECTOR ..........
  280    MP = 1
C
         DO 320 I = 1, UK
C
            DO 300 J = MP, UK
               RM1(I,J) = AR(I,J)
               RM2(I,J) = AI(I,J)
  300       CONTINUE
C
            RM1(I,I) = RM1(I,I) - RLAMBD
            RM2(I,I) = RM2(I,I) - ILAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C     .......... TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
         IF (UK .EQ. 1) GO TO 420
C
         DO 400 I = 2, UK
            MP = I - 1
            IF (PYTHAG(RM1(I,MP),RM2(I,MP)) .LE.
     1         PYTHAG(RM1(MP,MP),RM2(MP,MP))) GO TO 360
C
            DO 340 J = MP, UK
               Y = RM1(I,J)
               RM1(I,J) = RM1(MP,J)
               RM1(MP,J) = Y
               Y = RM2(I,J)
               RM2(I,J) = RM2(MP,J)
               RM2(MP,J) = Y
  340       CONTINUE
C
  360       IF (RM1(MP,MP) .EQ. 0.0E0 .AND. RM2(MP,MP) .EQ. 0.0E0)
     1         RM1(MP,MP) = EPS3
            CALL CDIV(RM1(I,MP),RM2(I,MP),RM1(MP,MP),RM2(MP,MP),X,Y)
            IF (X .EQ. 0.0E0 .AND. Y .EQ. 0.0E0) GO TO 400
C
            DO 380 J = I, UK
               RM1(I,J) = RM1(I,J) - X * RM1(MP,J) + Y * RM2(MP,J)
               RM2(I,J) = RM2(I,J) - X * RM2(MP,J) - Y * RM1(MP,J)
  380       CONTINUE
C
  400    CONTINUE
C
  420    IF (RM1(UK,UK) .EQ. 0.0E0 .AND. RM2(UK,UK) .EQ. 0.0E0)
     1      RM1(UK,UK) = EPS3
         ITS = 0
C     .......... BACK SUBSTITUTION
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0E0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C
            DO 680 J = IP1, UK
               X = X - RM1(I,J) * RV1(J) + RM2(I,J) * RV2(J)
               Y = Y - RM1(I,J) * RV2(J) - RM2(I,J) * RV1(J)
  680       CONTINUE
C
  700       CALL CDIV(X,Y,RM1(I,I),RM2(I,I),RV1(I),RV2(I))
  720    CONTINUE
C     .......... ACCEPTANCE TEST FOR EIGENVECTOR
C                AND NORMALIZATION ..........
         ITS = ITS + 1
         NORM = 0.0E0
         NORMV = 0.0E0
C
         DO 780 I = 1, UK
            X = PYTHAG(RV1(I),RV2(I))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
C
         IF (NORM .LT. GROWTO) GO TO 840
C     .......... ACCEPT VECTOR ..........
         X = RV1(J)
         Y = RV2(J)
C
         DO 820 I = 1, UK
            CALL CDIV(RV1(I),RV2(I),X,Y,ZR(I,S),ZI(I,S))
  820    CONTINUE
C
         IF (UK .EQ. N) GO TO 940
         J = UK + 1
         GO TO 900
C     .......... IN-LINE PROCEDURE FOR CHOOSING
C                A NEW STARTING VECTOR ..........
  840    IF (ITS .GE. UK) GO TO 880
         X = UKROOT
         Y = EPS3 / (X + 1.0E0)
         RV1(1) = EPS3
C
         DO 860 I = 2, UK
  860    RV1(I) = Y
C
         J = UK - ITS + 1
         RV1(J) = RV1(J) - EPS3 * X
         GO TO 660
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
         IERR = -K
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
            ZR(I,S) = 0.0E0
            ZI(I,S) = 0.0E0
  920    CONTINUE
C
  940    S = S + 1
  980 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1
      RETURN
      END
