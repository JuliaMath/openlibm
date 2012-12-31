*DECK BALANC
      SUBROUTINE BALANC (NM, N, A, LOW, IGH, SCALE)
C***BEGIN PROLOGUE  BALANC
C***PURPOSE  Balance a real general matrix and isolate eigenvalues
C            whenever possible.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1A
C***TYPE      SINGLE PRECISION (BALANC-S, CBAL-C)
C***KEYWORDS  EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure BALANCE,
C     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
C     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
C
C     This subroutine balances a REAL matrix and isolates
C     eigenvalues whenever possible.
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
C        A contains the input matrix to be balanced.  A is a
C          two-dimensional REAL array, dimensioned A(NM,N).
C
C     On OUTPUT
C
C        A contains the balanced matrix.
C
C        LOW and IGH are two INTEGER variables such that A(I,J)
C          is equal to zero if
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
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     The order in which the interchanges are made is N to IGH+1,
C     then 1 TO LOW-1.
C
C     Note that 1 is returned for IGH if IGH is zero formally.
C
C     The ALGOL procedure EXC contained in BALANCE appears in
C     BALANC  in line.  (Note that the ALGOL roles of identifiers
C     K,L have been reversed.)
C
C     Questions and comments should be directed to B. S. Garbow,
C     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
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
C***END PROLOGUE  BALANC
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,*),SCALE(*)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C***FIRST EXECUTABLE STATEMENT  BALANC
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
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
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
            IF (A(J,I) .NE. 0.0E0) GO TO 120
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
            IF (A(I,J) .NE. 0.0E0) GO TO 170
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
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
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
  250    A(I,J) = A(I,J) * G
C
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
