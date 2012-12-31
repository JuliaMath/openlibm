*DECK CBAL
      SUBROUTINE CBAL (NM, N, AR, AI, LOW, IGH, SCALE)
C***BEGIN PROLOGUE  CBAL
C***PURPOSE  Balance a complex general matrix and isolate eigenvalues
C            whenever possible.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1A
C***TYPE      COMPLEX (BALANC-S, CBAL-C)
C***KEYWORDS  EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure
C     CBALANCE, which is a complex version of BALANCE,
C     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     This subroutine balances a COMPLEX matrix and isolates
C     eigenvalues whenever possible.
C
C     On INPUT
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, AR and AI, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A=(AR,AI).  N is an INTEGER
C          variable.  N must be less than or equal to NM.
C
C        AR and AI contain the real and imaginary parts,
C          respectively, of the complex matrix to be balanced.
C          AR and AI are two-dimensional REAL arrays, dimensioned
C          AR(NM,N) and AI(NM,N).
C
C     On OUTPUT
C
C        AR and AI contain the real and imaginary parts,
C          respectively, of the balanced matrix.
C
C        LOW and IGH are two INTEGER variables such that AR(I,J)
C          and AI(I,J) are equal to zero if
C           (1) I is greater than J and
C           (2) J=1,...,LOW-1 or I=IGH+1,...,N.
C
C        SCALE contains information determining the permutations and
C          scaling factors used.  SCALE is a one-dimensional REAL array,
C          dimensioned SCALE(N).
C
C     Suppose that the principal submatrix in rows LOW through IGH
C     has been balanced, that P(J) denotes the index interchanged
C     with J during the permutation step, and that the elements
C     of the diagonal matrix used are denoted by D(I,J).  Then
C        SCALE(J) = P(J),    for J = 1,...,LOW-1
C                 = D(J,J)       J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     The order in which the interchanges are made is N to IGH+1,
C     then 1 to LOW-1.
C
C     Note that 1 is returned for IGH if IGH is zero formally.
C
C     The ALGOL procedure EXC contained in CBALANCE appears in
C     CBAL  in line.  (Note that the ALGOL roles of identifiers
C     K,L have been reversed.)
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CBAL
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL AR(NM,*),AI(NM,*),SCALE(*)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C     THE FOLLOWING PORTABLE VALUE OF RADIX WORKS WELL ENOUGH
C     FOR ALL MACHINES WHOSE BASE IS A POWER OF TWO.
C
C***FIRST EXECUTABLE STATEMENT  CBAL
      RADIX = 16
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     .......... IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (AR(J,I) .NE. 0.0E0 .OR. AI(J,I) .NE. 0.0E0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ..........
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (AR(I,J) .NE. 0.0E0 .OR. AI(I,J) .NE. 0.0E0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0E0
C     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0E0
         R = 0.0E0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))
  200    CONTINUE
C     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
         IF (C .EQ. 0.0E0 .OR. R .EQ. 0.0E0) GO TO 270
         G = R / RADIX
         F = 1.0E0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     .......... NOW BALANCE ..........
  240    IF ((C + R) / F .GE. 0.95E0 * S) GO TO 270
         G = 1.0E0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
C
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
