*DECK TSTURM
      SUBROUTINE TSTURM (NM, N, EPS1, D, E, E2, LB, UB, MM, M, W, Z,
     +   IERR, RV1, RV2, RV3, RV4, RV5, RV6)
C***BEGIN PROLOGUE  TSTURM
C***PURPOSE  Find those eigenvalues of a symmetric tridiagonal matrix
C            in a given interval and their associated eigenvectors by
C            Sturm sequencing.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (TSTURM-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine finds those eigenvalues of a TRIDIAGONAL
C     SYMMETRIC matrix which lie in a specified interval and their
C     associated eigenvectors, using bisection and inverse iteration.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        EPS1 is an absolute error tolerance for the computed eigen-
C          values.  It should be chosen so that the accuracy of these
C          eigenvalues is commensurate with relative perturbations of
C          the order of the relative machine precision in the matrix
C          elements.  If the input EPS1 is non-positive, it is reset
C          for each submatrix to a default value, namely, minus the
C          product of the relative machine precision and the 1-norm of
C          the submatrix.  EPS1 is a REAL variable.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional REAL array, dimensioned
C          E(N).
C
C        E2 contains the squares of the corresponding elements of E.
C          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
C          dimensioned E2(N).
C
C        LB and UB define the interval to be searched for eigenvalues.
C          If LB is not less than UB, no eigenvalues will be found.
C          LB and UB are REAL variables.
C
C        MM should be set to an upper bound for the number of
C          eigenvalues in the interval.  MM is an INTEGER variable.
C          WARNING -  If more than MM eigenvalues are determined to lie
C          in the interval, an error return is made with no values or
C          vectors found.
C
C     On Output
C
C        EPS1 is unaltered unless it has been reset to its
C          (last) default value.
C
C        D and E are unaltered.
C
C        Elements of E2, corresponding to elements of E regarded as
C          negligible, have been replaced by zero causing the matrix to
C          split into a direct sum of submatrices.  E2(1) is also set
C          to zero.
C
C        M is the number of eigenvalues determined to lie in (LB,UB).
C          M is an INTEGER variable.
C
C        W contains the M eigenvalues in ascending order if the matrix
C          does not split.  If the matrix splits, the eigenvalues are
C          in ascending order for each submatrix.  If a vector error
C          exit is made, W contains those values already found.  W is a
C          one-dimensional REAL array, dimensioned W(MM).
C
C        Z contains the associated set of orthonormal eigenvectors.
C          If an error exit is made, Z contains those vectors already
C          found.  Z is a one-dimensional REAL array, dimensioned
C          Z(NM,MM).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          3*N+1      if M exceeds MM no eigenvalues or eigenvectors
C                     are computed,
C          4*N+J      if the eigenvector corresponding to the J-th
C                     eigenvalue fails to converge in 5 iterations, then
C                     the eigenvalues and eigenvectors in W and Z should
C                     be correct for indices 1, 2, ..., J-1.
C
C        RV1, RV2, RV3, RV4, RV5, and RV6 are temporary storage arrays,
C          dimensioned RV1(N), RV2(N), RV3(N), RV4(N), RV5(N), and
C          RV6(N).
C
C     The ALGOL procedure STURMCNT contained in TRISTURM
C     appears in TSTURM in-line.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  TSTURM
C
      INTEGER I,J,K,M,N,P,Q,R,S,II,IP,JJ,MM,M1,M2,NM,ITS
      INTEGER IERR,GROUP,ISTURM
      REAL D(*),E(*),E2(*),W(*),Z(NM,*)
      REAL RV1(*),RV2(*),RV3(*),RV4(*),RV5(*),RV6(*)
      REAL U,V,LB,T1,T2,UB,UK,XU,X0,X1,EPS1,EPS2,EPS3,EPS4
      REAL NORM,MACHEP,S1,S2
      LOGICAL FIRST
C
      SAVE FIRST, MACHEP
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  TSTURM
      IF (FIRST) THEN
         MACHEP = R1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      IERR = 0
      T1 = LB
      T2 = UB
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO 40 I = 1, N
         IF (I .EQ. 1) GO TO 20
         S1 = ABS(D(I)) + ABS(D(I-1))
         S2 = S1 + ABS(E(I))
         IF (S2 .GT. S1) GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
C     .......... DETERMINE THE NUMBER OF EIGENVALUES
C                IN THE INTERVAL ..........
      P = 1
      Q = N
      X1 = UB
      ISTURM = 1
      GO TO 320
   60 M = S
      X1 = LB
      ISTURM = 2
      GO TO 320
   80 M = M - S
      IF (M .GT. MM) GO TO 980
      Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0E0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0E0
         V = 0.0E0
         IF (Q .EQ. N) GO TO 110
         U = ABS(E(Q+1))
         V = E2(Q+1)
  110    XU = MIN(D(Q)-(X1+U),XU)
         X0 = MAX(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0E0) GO TO 140
  120 CONTINUE
C
  140 X1 = MAX(ABS(XU),ABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0E0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      R = R + 1
C
      DO 160 I = 1, N
  160 Z(I,R) = 0.0E0
C
      W(R) = D(P)
      Z(P,R) = 1.0E0
      GO TO 940
  180 X1 = X1 * (Q-P+1)
      LB = MAX(T1,XU-X1)
      UB = MIN(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     .......... FIND ROOTS BY BISECTION ..........
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     .......... LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
      K = M2
  250    XU = LB
C     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     .......... NEXT BISECTION STEP ..........
  300    X1 = (XU + X0) * 0.5E0
         S1 = 2.0E0*(ABS(XU) + ABS(X0) + ABS(EPS1))
         S2 = S1 + ABS(X0 - XU)
         IF (S2 .EQ. S1) GO TO 420
C     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
  320    S = P - 1
         U = 1.0E0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0E0) GO TO 325
            V = ABS(E(I)) / MACHEP
            IF (E2(I) .EQ. 0.0E0) V = 0.0E0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0E0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     .......... REFINE INTERVALS ..........
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     .......... K-TH EIGENVALUE FOUND ..........
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     .......... FIND VECTORS BY INVERSE ITERATION ..........
      NORM = ABS(D(P))
      IP = P + 1
C
      DO 500 I = IP, Q
  500 NORM = MAX(NORM, ABS(D(I)) + ABS(E(I)))
C     .......... EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
      EPS2 = 1.0E-3 * NORM
      UK = SQRT(REAL(Q-P+5))
      EPS3 = UK * MACHEP * NORM
      EPS4 = UK * EPS3
      UK = EPS4 / SQRT(UK)
      GROUP = 0
      S = P
C
      DO 920 K = M1, M2
         R = R + 1
         ITS = 1
         W(R) = RV5(K)
         X1 = RV5(K)
C     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
         IF (K .EQ. M1) GO TO 520
         IF (X1 - X0 .GE. EPS2) GROUP = -1
         GROUP = GROUP + 1
         IF (X1 .LE. X0) X1 = X0 + EPS3
C     .......... ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (ABS(E(I)) .LT. ABS(U)) GO TO 540
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0E0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. 0.0E0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0E0
         RV3(Q) = 0.0E0
C     .......... BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ..........
         IF (GROUP .EQ. 0) GO TO 700
C
         DO 680 JJ = 1, GROUP
            J = R - GROUP - 1 + JJ
            XU = 0.0E0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0E0
C
         DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
C
         IF (NORM .GE. 1.0E0) GO TO 840
C     .......... FORWARD SUBSTITUTION ..........
         IF (ITS .EQ. 5) GO TO 960
         IF (NORM .NE. 0.0E0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ..........
  780    DO 820 I = IP, Q
            U = RV6(I)
C     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ..........
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     .......... NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
C
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
C
         XU = 1.0E0 / SQRT(U)
C
         DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  960 IERR = 4 * N + R
      GO TO 1001
C     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
C                EIGENVALUES IN INTERVAL ..........
  980 IERR = 3 * N + 1
 1001 LB = T1
      UB = T2
      RETURN
      END
