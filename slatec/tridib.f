*DECK TRIDIB
      SUBROUTINE TRIDIB (N, EPS1, D, E, E2, LB, UB, M11, M, W, IND,
     +   IERR, RV4, RV5)
C***BEGIN PROLOGUE  TRIDIB
C***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix
C            in a given interval using Sturm sequencing.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (TRIDIB-S)
C***KEYWORDS  EIGENVALUES OF A REAL SYMMETRIC MATRIX, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BISECT,
C     NUM. MATH. 9, 386-393(1967) by Barth, Martin, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     This subroutine finds those eigenvalues of a TRIDIAGONAL
C     SYMMETRIC matrix between specified boundary indices,
C     using bisection.
C
C     On Input
C
C        N is the order of the matrix.  N is an INTEGER variable.
C
C        EPS1 is an absolute error tolerance for the computed eigen-
C          values.  If the input EPS1 is non-positive, it is reset for
C          each submatrix to a default value, namely, minus the product
C          of the relative machine precision and the 1-norm of the
C          submatrix.  EPS1 is a REAL variable.
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
C        M11 specifies the lower boundary index for the set of desired
C          eigenvalues.  M11 is an INTEGER variable.
C
C        M specifies the number of eigenvalues desired.  The upper
C          boundary index M22 is then obtained as M22=M11+M-1.
C          M is an INTEGER variable.
C
C     On Output
C
C        EPS1 is unaltered unless it has been reset to its
C          (last) default value.
C
C        D and E are unaltered.
C
C        Elements of E2, corresponding to elements of E regarded
C          as negligible, have been replaced by zero causing the
C          matrix to split into a direct sum of submatrices.
C          E2(1) is also set to zero.
C
C        LB and UB define an interval containing exactly the desired
C          eigenvalues.  LB and UB are REAL variables.
C
C        W contains, in its first M positions, the eigenvalues
C          between indices M11 and M22 in ascending order.
C          W is a one-dimensional REAL array, dimensioned W(M).
C
C        IND contains in its first M positions the submatrix indices
C          associated with the corresponding eigenvalues in W --
C          1 for eigenvalues belonging to the first submatrix from
C          the top, 2 for those belonging to the second submatrix, etc.
C          IND is an one-dimensional INTEGER array, dimensioned IND(M).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          3*N+1      if multiple eigenvalues at index M11 make
C                     unique selection of LB impossible,
C          3*N+2      if multiple eigenvalues at index M22 make
C                     unique selection of UB impossible.
C
C        RV4 and RV5 are one-dimensional REAL arrays used for temporary
C          storage of the lower and upper bounds for the eigenvalues in
C          the bisection process.  RV4 and RV5 are dimensioned RV4(N)
C          and RV5(N).
C
C     Note that subroutine TQL1, IMTQL1, or TQLRAT is generally faster
C     than TRIDIB, if more than N/4 eigenvalues are to be found.
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
C***END PROLOGUE  TRIDIB
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      REAL D(*),E(*),E2(*),W(*),RV4(*),RV5(*)
      REAL U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP,S1,S2
      INTEGER IND(*)
      LOGICAL FIRST
C
      SAVE FIRST, MACHEP
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  TRIDIB
      IF (FIRST) THEN
         MACHEP = R1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
      DO 40 I = 1, N
         X1 = U
         U = 0.0E0
         IF (I .NE. N) U = ABS(E(I+1))
         XU = MIN(D(I)-(X1+U),XU)
         X0 = MAX(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         S1 = ABS(D(I)) + ABS(D(I-1))
         S2 = S1 + ABS(E(I))
         IF (S2 .GT. S1) GO TO 40
   20    E2(I) = 0.0E0
   40 CONTINUE
C
      X1 = MAX(ABS(XU),ABS(X0)) * MACHEP * N
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     .......... DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES ..........
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5E0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
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
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
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
         S1 = ABS(XU) + ABS(X0) + ABS(EPS1)
         S2 = S1 + ABS(X0-XU)/2.0E0
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
C     .......... ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ..........
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     .......... SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES ..........
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
      END
