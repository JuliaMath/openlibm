*DECK RATQR
      SUBROUTINE RATQR (N, EPS1, D, E, E2, M, W, IND, BD, TYPE, IDEF,
     +   IERR)
C***BEGIN PROLOGUE  RATQR
C***PURPOSE  Compute the largest or smallest eigenvalues of a symmetric
C            tridiagonal matrix using the rational QR method with Newton
C            correction.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (RATQR-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure RATQR,
C     NUM. MATH. 11, 264-272(1968) by REINSCH and BAUER.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 257-265(1971).
C
C     This subroutine finds the algebraically smallest or largest
C     eigenvalues of a SYMMETRIC TRIDIAGONAL matrix by the
C     rational QR method with Newton corrections.
C
C     On Input
C
C        N is the order of the matrix.  N is an INTEGER variable.
C
C        EPS1 is a theoretical absolute error tolerance for the
C          computed eigenvalues.  If the input EPS1 is non-positive, or
C          indeed smaller than its default value, it is reset at each
C          iteration to the respective default value, namely, the
C          product of the relative machine precision and the magnitude
C          of the current eigenvalue iterate.  The theoretical absolute
C          error in the K-th eigenvalue is usually not greater than
C          K times EPS1.  EPS1 is a REAL variable.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional REAL array, dimensioned
C          E(N).
C
C        E2 contains the squares of the corresponding elements of E in
C          its last N-1 positions.  E2(1) is arbitrary.  E2 is a one-
C          dimensional REAL array, dimensioned E2(N).
C
C        M is the number of eigenvalues to be found.  M is an INTEGER
C          variable.
C
C        IDEF should be set to 1 if the input matrix is known to be
C          positive definite, to -1 if the input matrix is known to
C          be negative definite, and to 0 otherwise.  IDEF is an
C          INTEGER variable.
C
C        TYPE should be set to .TRUE. if the smallest eigenvalues are
C          to be found, and to .FALSE. if the largest eigenvalues are
C          to be found.  TYPE is a LOGICAL variable.
C
C     On Output
C
C        EPS1 is unaltered unless it has been reset to its
C          (last) default value.
C
C        D and E are unaltered (unless W overwrites D).
C
C        Elements of E2, corresponding to elements of E regarded as
C          negligible, have been replaced by zero causing the matrix
C          to split into a direct sum of submatrices.  E2(1) is set
C          to 0.0e0 if the smallest eigenvalues have been found, and
C          to 2.0e0 if the largest eigenvalues have been found.  E2
C          is otherwise unaltered (unless overwritten by BD).
C
C        W contains the M algebraically smallest eigenvalues in
C          ascending order, or the M largest eigenvalues in descending
C          order.  If an error exit is made because of an incorrect
C          specification of IDEF, no eigenvalues are found.  If the
C          Newton iterates for a particular eigenvalue are not monotone,
C          the best estimate obtained is returned and IERR is set.
C          W is a one-dimensional REAL array, dimensioned W(N).  W need
C          not be distinct from D.
C
C        IND contains in its first M positions the submatrix indices
C          associated with the corresponding eigenvalues in W --
C          1 for eigenvalues belonging to the first submatrix from
C          the top, 2 for those belonging to the second submatrix, etc.
C          IND is an one-dimensional INTEGER array, dimensioned IND(N).
C
C        BD contains refined bounds for the theoretical errors of the
C          corresponding eigenvalues in W.  These bounds are usually
C          within the tolerance specified by EPS1.  BD is a one-
C          dimensional REAL array, dimensioned BD(N).  BD need not be
C          distinct from E2.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          6*N+1      if  IDEF  is set to 1 and  TYPE  to .TRUE.
C                     when the matrix is NOT positive definite, or
C                     if  IDEF  is set to -1 and  TYPE  to .FALSE.
C                     when the matrix is NOT negative definite,
C                     no eigenvalues are computed, or
C                     M is greater than N,
C          5*N+K      if successive iterates to the K-th eigenvalue
C                     are NOT monotone increasing, where K refers
C                     to the last such occurrence.
C
C     Note that subroutine TRIDIB is generally faster and more
C     accurate than RATQR if the eigenvalues are clustered.
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
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RATQR
C
      INTEGER I,J,K,M,N,II,JJ,K1,IDEF,IERR,JDEF
      REAL D(*),E(*),E2(*),W(*),BD(*)
      REAL F,P,Q,R,S,EP,QP,ERR,TOT,EPS1,DELTA,MACHEP
      INTEGER IND(*)
      LOGICAL FIRST, TYPE
C
      SAVE FIRST, MACHEP
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  RATQR
      IF (FIRST) THEN
         MACHEP = R1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      IERR = 0
      JDEF = IDEF
C     .......... COPY D ARRAY INTO W ..........
      DO 20 I = 1, N
   20 W(I) = D(I)
C
      IF (TYPE) GO TO 40
      J = 1
      GO TO 400
   40 ERR = 0.0E0
      S = 0.0E0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DEFINE
C                INITIAL SHIFT FROM LOWER GERSCHGORIN BOUND.
C                COPY E2 ARRAY INTO BD ..........
      TOT = W(1)
      Q = 0.0E0
      J = 0
C
      DO 100 I = 1, N
         P = Q
         IF (I .EQ. 1) GO TO 60
         IF (P .GT. MACHEP * (ABS(D(I)) + ABS(D(I-1)))) GO TO 80
   60    E2(I) = 0.0E0
   80    BD(I) = E2(I)
C     .......... COUNT ALSO IF ELEMENT OF E2 HAS UNDERFLOWED ..........
         IF (E2(I) .EQ. 0.0E0) J = J + 1
         IND(I) = J
         Q = 0.0E0
         IF (I .NE. N) Q = ABS(E(I+1))
         TOT = MIN(W(I)-P-Q,TOT)
  100 CONTINUE
C
      IF (JDEF .EQ. 1 .AND. TOT .LT. 0.0E0) GO TO 140
C
      DO 110 I = 1, N
  110 W(I) = W(I) - TOT
C
      GO TO 160
  140 TOT = 0.0E0
C
  160 DO 360 K = 1, M
C     .......... NEXT QR TRANSFORMATION ..........
  180    TOT = TOT + S
         DELTA = W(N) - S
         I = N
         F = ABS(MACHEP*TOT)
         IF (EPS1 .LT. F) EPS1 = F
         IF (DELTA .GT. EPS1) GO TO 190
         IF (DELTA .LT. (-EPS1)) GO TO 1000
         GO TO 300
C     .......... REPLACE SMALL SUB-DIAGONAL SQUARES BY ZERO
C                TO REDUCE THE INCIDENCE OF UNDERFLOWS ..........
  190    IF (K .EQ. N) GO TO 210
         K1 = K + 1
         DO 200 J = K1, N
            IF (BD(J) .LE. (MACHEP*(W(J)+W(J-1))) ** 2) BD(J) = 0.0E0
  200    CONTINUE
C
  210    F = BD(N) / DELTA
         QP = DELTA + F
         P = 1.0E0
         IF (K .EQ. N) GO TO 260
         K1 = N - K
C     .......... FOR I=N-1 STEP -1 UNTIL K DO -- ..........
         DO 240 II = 1, K1
            I = N - II
            Q = W(I) - S - F
            R = Q / QP
            P = P * R + 1.0E0
            EP = F * R
            W(I+1) = QP + EP
            DELTA = Q - EP
            IF (DELTA .GT. EPS1) GO TO 220
            IF (DELTA .LT. (-EPS1)) GO TO 1000
            GO TO 300
  220       F = BD(I) / Q
            QP = DELTA + F
            BD(I+1) = QP * EP
  240    CONTINUE
C
  260    W(K) = QP
         S = QP / P
         IF (TOT + S .GT. TOT) GO TO 180
C     .......... SET ERROR -- IRREGULAR END OF ITERATION.
C                DEFLATE MINIMUM DIAGONAL ELEMENT ..........
         IERR = 5 * N + K
         S = 0.0E0
         DELTA = QP
C
         DO 280 J = K, N
            IF (W(J) .GT. DELTA) GO TO 280
            I = J
            DELTA = W(J)
  280    CONTINUE
C     .......... CONVERGENCE ..........
  300    IF (I .LT. N) BD(I+1) = BD(I) * F / QP
         II = IND(I)
         IF (I .EQ. K) GO TO 340
         K1 = I - K
C     .......... FOR J=I-1 STEP -1 UNTIL K DO -- ..........
         DO 320 JJ = 1, K1
            J = I - JJ
            W(J+1) = W(J) - S
            BD(J+1) = BD(J)
            IND(J+1) = IND(J)
  320    CONTINUE
C
  340    W(K) = TOT
         ERR = ERR + ABS(DELTA)
         BD(K) = ERR
         IND(K) = II
  360 CONTINUE
C
      IF (TYPE) GO TO 1001
      F = BD(1)
      E2(1) = 2.0E0
      BD(1) = F
      J = 2
C     .......... NEGATE ELEMENTS OF W FOR LARGEST VALUES ..........
  400 DO 500 I = 1, N
  500 W(I) = -W(I)
C
      JDEF = -JDEF
      GO TO (40,1001), J
C     .......... SET ERROR -- IDEF SPECIFIED INCORRECTLY ..........
 1000 IERR = 6 * N + 1
 1001 RETURN
      END
