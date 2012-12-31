*DECK INVIT
      SUBROUTINE INVIT (NM, N, A, WR, WI, SELECT, MM, M, Z, IERR, RM1,
     +   RV1, RV2)
C***BEGIN PROLOGUE  INVIT
C***PURPOSE  Compute the eigenvectors of a real upper Hessenberg
C            matrix associated with specified eigenvalues by inverse
C            iteration.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C2B
C***TYPE      SINGLE PRECISION (INVIT-S, CINVIT-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure INVIT
C     by Peters and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     This subroutine finds those eigenvectors of a REAL UPPER
C     Hessenberg matrix corresponding to specified eigenvalues,
C     using inverse iteration.
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
C        A contains the upper Hessenberg matrix.  A is a two-dimensional
C          REAL array, dimensioned A(NM,N).
C
C        WR and WI contain the real and imaginary parts, respectively,
C          of the eigenvalues of the Hessenberg matrix.  The eigenvalues
C          must be stored in a manner identical to that output by
C          subroutine  HQR,  which recognizes possible splitting of the
C          matrix.  WR and WI are one-dimensional REAL arrays,
C          dimensioned WR(N) and WI(N).
C
C        SELECT specifies the eigenvectors to be found. The
C          eigenvector corresponding to the J-th eigenvalue is
C          specified by setting SELECT(J) to .TRUE.  SELECT is a
C          one-dimensional LOGICAL array, dimensioned SELECT(N).
C
C        MM should be set to an upper bound for the number of
C          columns required to store the eigenvectors to be found.
C          NOTE that two columns are required to store the
C          eigenvector corresponding to a complex eigenvalue.  One
C          column is required to store the eigenvector corresponding
C          to a real eigenvalue.  MM is an INTEGER variable.
C
C     On OUTPUT
C
C        A and WI are unaltered.
C
C        WR may have been altered since close eigenvalues are perturbed
C          slightly in searching for independent eigenvectors.
C
C        SELECT may have been altered.  If the elements corresponding
C          to a pair of conjugate complex eigenvalues were each
C          initially set to .TRUE., the program resets the second of
C          the two elements to .FALSE.
C
C        M is the number of columns actually used to store the
C          eigenvectors.  M is an INTEGER variable.
C
C        Z contains the real and imaginary parts of the eigenvectors.
C          The eigenvectors are packed into the columns of Z starting
C          at the first column.  If the next selected eigenvalue is
C          real, the next column of Z contains its eigenvector.  If the
C          eigenvalue is complex, the next two columns of Z contain the
C          real and imaginary parts of its eigenvector, with the real
C          part first.  The eigenvectors are normalized so that the
C          component of largest magnitude is 1. Any vector which fails
C          the acceptance test is set to zero.  Z is a two-dimensional
C          REAL array, dimensioned Z(NM,MM).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          -(2*N+1)   if more than MM columns of Z are necessary
C                     to store the eigenvectors corresponding to
C                     the specified eigenvalues (in this case, M is
C                     equal to the number of columns of Z containing
C                     eigenvectors already computed),
C          -K         if the iteration corresponding to the K-th
C                     value fails (if this occurs more than once, K
C                     is the index of the last occurrence); the
C                     corresponding columns of Z are set to zero
C                     vectors,
C          -(N+K)     if both error situations occur.
C
C        RM1 is a two-dimensional REAL array used for temporary storage.
C          This array holds the triangularized form of the upper
C          Hessenberg matrix used in the inverse iteration process.
C          RM1 is dimensioned RM1(N,N).
C
C        RV1 and RV2 are one-dimensional REAL arrays used for temporary
C          storage.  They hold the approximate eigenvectors during the
C          inverse iteration process.  RV1 and RV2 are dimensioned
C          RV1(N) and RV2(N).
C
C     The ALGOL procedure GUESSVEC appears in INVIT in-line.
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
C***END PROLOGUE  INVIT
C
      INTEGER I,J,K,L,M,N,S,II,IP,MM,MP,NM,NS,N1,UK,IP1,ITS,KM1,IERR
      REAL A(NM,*),WR(*),WI(*),Z(NM,*)
      REAL RM1(N,*),RV1(*),RV2(*)
      REAL T,W,X,Y,EPS3
      REAL NORM,NORMV,GROWTO,ILAMBD,RLAMBD,UKROOT
      REAL PYTHAG
      LOGICAL SELECT(N)
C
C***FIRST EXECUTABLE STATEMENT  INVIT
      IERR = 0
      UK = 0
      S = 1
C     .......... IP = 0, REAL EIGENVALUE
C                     1, FIRST OF CONJUGATE COMPLEX PAIR
C                    -1, SECOND OF CONJUGATE COMPLEX PAIR ..........
      IP = 0
      N1 = N - 1
C
      DO 980 K = 1, N
         IF (WI(K) .EQ. 0.0E0 .OR. IP .LT. 0) GO TO 100
         IP = 1
         IF (SELECT(K) .AND. SELECT(K+1)) SELECT(K+1) = .FALSE.
  100    IF (.NOT. SELECT(K)) GO TO 960
         IF (WI(K) .NE. 0.0E0) S = S + 1
         IF (S .GT. MM) GO TO 1000
         IF (UK .GE. K) GO TO 200
C     .......... CHECK FOR POSSIBLE SPLITTING ..........
         DO 120 UK = K, N
            IF (UK .EQ. N) GO TO 140
            IF (A(UK+1,UK) .EQ. 0.0E0) GO TO 140
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
  160       X = X + ABS(A(I,J))
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
C     .......... GROWTO IS THE CRITERION FOR THE GROWTH ..........
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
C     .......... PERTURB CONJUGATE EIGENVALUE TO MATCH ..........
         IP1 = K + IP
         WR(IP1) = RLAMBD
C     .......... FORM UPPER HESSENBERG A-RLAMBD*I (TRANSPOSED)
C                AND INITIAL REAL VECTOR ..........
  280    MP = 1
C
         DO 320 I = 1, UK
C
            DO 300 J = MP, UK
  300       RM1(J,I) = A(I,J)
C
            RM1(I,I) = RM1(I,I) - RLAMBD
            MP = I
            RV1(I) = EPS3
  320    CONTINUE
C
         ITS = 0
         IF (ILAMBD .NE. 0.0E0) GO TO 520
C     .......... REAL EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3 ..........
         IF (UK .EQ. 1) GO TO 420
C
         DO 400 I = 2, UK
            MP = I - 1
            IF (ABS(RM1(MP,I)) .LE. ABS(RM1(MP,MP))) GO TO 360
C
            DO 340 J = MP, UK
               Y = RM1(J,I)
               RM1(J,I) = RM1(J,MP)
               RM1(J,MP) = Y
  340       CONTINUE
C
  360       IF (RM1(MP,MP) .EQ. 0.0E0) RM1(MP,MP) = EPS3
            X = RM1(MP,I) / RM1(MP,MP)
            IF (X .EQ. 0.0E0) GO TO 400
C
            DO 380 J = I, UK
  380       RM1(J,I) = RM1(J,I) - X * RM1(J,MP)
C
  400    CONTINUE
C
  420    IF (RM1(UK,UK) .EQ. 0.0E0) RM1(UK,UK) = EPS3
C     .......... BACK SUBSTITUTION FOR REAL VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  440    DO 500 II = 1, UK
            I = UK + 1 - II
            Y = RV1(I)
            IF (I .EQ. UK) GO TO 480
            IP1 = I + 1
C
            DO 460 J = IP1, UK
  460       Y = Y - RM1(J,I) * RV1(J)
C
  480       RV1(I) = Y / RM1(I,I)
  500    CONTINUE
C
         GO TO 740
C     .......... COMPLEX EIGENVALUE.
C                TRIANGULAR DECOMPOSITION WITH INTERCHANGES,
C                REPLACING ZERO PIVOTS BY EPS3.  STORE IMAGINARY
C                PARTS IN UPPER TRIANGLE STARTING AT (1,3) ..........
  520    NS = N - S
         Z(1,S-1) = -ILAMBD
         Z(1,S) = 0.0E0
         IF (N .EQ. 2) GO TO 550
         RM1(1,3) = -ILAMBD
         Z(1,S-1) = 0.0E0
         IF (N .EQ. 3) GO TO 550
C
         DO 540 I = 4, N
  540    RM1(1,I) = 0.0E0
C
  550    DO 640 I = 2, UK
            MP = I - 1
            W = RM1(MP,I)
            IF (I .LT. N) T = RM1(MP,I+1)
            IF (I .EQ. N) T = Z(MP,S-1)
            X = RM1(MP,MP) * RM1(MP,MP) + T * T
            IF (W * W .LE. X) GO TO 580
            X = RM1(MP,MP) / W
            Y = T / W
            RM1(MP,MP) = W
            IF (I .LT. N) RM1(MP,I+1) = 0.0E0
            IF (I .EQ. N) Z(MP,S-1) = 0.0E0
C
            DO 560 J = I, UK
               W = RM1(J,I)
               RM1(J,I) = RM1(J,MP) - X * W
               RM1(J,MP) = W
               IF (J .LT. N1) GO TO 555
               L = J - NS
               Z(I,L) = Z(MP,L) - Y * W
               Z(MP,L) = 0.0E0
               GO TO 560
  555          RM1(I,J+2) = RM1(MP,J+2) - Y * W
               RM1(MP,J+2) = 0.0E0
  560       CONTINUE
C
            RM1(I,I) = RM1(I,I) - Y * ILAMBD
            IF (I .LT. N1) GO TO 570
            L = I - NS
            Z(MP,L) = -ILAMBD
            Z(I,L) = Z(I,L) + X * ILAMBD
            GO TO 640
  570       RM1(MP,I+2) = -ILAMBD
            RM1(I,I+2) = RM1(I,I+2) + X * ILAMBD
            GO TO 640
  580       IF (X .NE. 0.0E0) GO TO 600
            RM1(MP,MP) = EPS3
            IF (I .LT. N) RM1(MP,I+1) = 0.0E0
            IF (I .EQ. N) Z(MP,S-1) = 0.0E0
            T = 0.0E0
            X = EPS3 * EPS3
  600       W = W / X
            X = RM1(MP,MP) * W
            Y = -T * W
C
            DO 620 J = I, UK
               IF (J .LT. N1) GO TO 610
               L = J - NS
               T = Z(MP,L)
               Z(I,L) = -X * T - Y * RM1(J,MP)
               GO TO 615
  610          T = RM1(MP,J+2)
               RM1(I,J+2) = -X * T - Y * RM1(J,MP)
  615          RM1(J,I) = RM1(J,I) - X * RM1(J,MP) + Y * T
  620       CONTINUE
C
            IF (I .LT. N1) GO TO 630
            L = I - NS
            Z(I,L) = Z(I,L) - ILAMBD
            GO TO 640
  630       RM1(I,I+2) = RM1(I,I+2) - ILAMBD
  640    CONTINUE
C
         IF (UK .LT. N1) GO TO 650
         L = UK - NS
         T = Z(UK,L)
         GO TO 655
  650    T = RM1(UK,UK+2)
  655    IF (RM1(UK,UK) .EQ. 0.0E0 .AND. T .EQ. 0.0E0) RM1(UK,UK) = EPS3
C     .......... BACK SUBSTITUTION FOR COMPLEX VECTOR
C                FOR I=UK STEP -1 UNTIL 1 DO -- ..........
  660    DO 720 II = 1, UK
            I = UK + 1 - II
            X = RV1(I)
            Y = 0.0E0
            IF (I .EQ. UK) GO TO 700
            IP1 = I + 1
C
            DO 680 J = IP1, UK
               IF (J .LT. N1) GO TO 670
               L = J - NS
               T = Z(I,L)
               GO TO 675
  670          T = RM1(I,J+2)
  675          X = X - RM1(J,I) * RV1(J) + T * RV2(J)
               Y = Y - RM1(J,I) * RV2(J) - T * RV1(J)
  680       CONTINUE
C
  700       IF (I .LT. N1) GO TO 710
            L = I - NS
            T = Z(I,L)
            GO TO 715
  710       T = RM1(I,I+2)
  715       CALL CDIV(X,Y,RM1(I,I),T,RV1(I),RV2(I))
  720    CONTINUE
C     .......... ACCEPTANCE TEST FOR REAL OR COMPLEX
C                EIGENVECTOR AND NORMALIZATION ..........
  740    ITS = ITS + 1
         NORM = 0.0E0
         NORMV = 0.0E0
C
         DO 780 I = 1, UK
            IF (ILAMBD .EQ. 0.0E0) X = ABS(RV1(I))
            IF (ILAMBD .NE. 0.0E0) X = PYTHAG(RV1(I),RV2(I))
            IF (NORMV .GE. X) GO TO 760
            NORMV = X
            J = I
  760       NORM = NORM + X
  780    CONTINUE
C
         IF (NORM .LT. GROWTO) GO TO 840
C     .......... ACCEPT VECTOR ..........
         X = RV1(J)
         IF (ILAMBD .EQ. 0.0E0) X = 1.0E0 / X
         IF (ILAMBD .NE. 0.0E0) Y = RV2(J)
C
         DO 820 I = 1, UK
            IF (ILAMBD .NE. 0.0E0) GO TO 800
            Z(I,S) = RV1(I) * X
            GO TO 820
  800       CALL CDIV(RV1(I),RV2(I),X,Y,Z(I,S-1),Z(I,S))
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
         IF (ILAMBD .EQ. 0.0E0) GO TO 440
         GO TO 660
C     .......... SET ERROR -- UNACCEPTED EIGENVECTOR ..........
  880    J = 1
         IERR = -K
C     .......... SET REMAINING VECTOR COMPONENTS TO ZERO ..........
  900    DO 920 I = J, N
            Z(I,S) = 0.0E0
            IF (ILAMBD .NE. 0.0E0) Z(I,S-1) = 0.0E0
  920    CONTINUE
C
  940    S = S + 1
  960    IF (IP .EQ. (-1)) IP = 0
         IF (IP .EQ. 1) IP = -1
  980 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF EIGENVECTOR
C                SPACE REQUIRED ..........
 1000 IF (IERR .NE. 0) IERR = IERR - N
      IF (IERR .EQ. 0) IERR = -(2 * N + 1)
 1001 M = S - 1 - ABS(IP)
      RETURN
      END
