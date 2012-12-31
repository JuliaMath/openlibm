*DECK DHFTI
      SUBROUTINE DHFTI (A, MDA, M, N, B, MDB, NB, TAU, KRANK, RNORM, H,
     +   G, IP)
C***BEGIN PROLOGUE  DHFTI
C***PURPOSE  Solve a least squares problem for banded matrices using
C            sequential accumulation of rows of the data matrix.
C            Exactly one right-hand side vector is permitted.
C***LIBRARY   SLATEC
C***CATEGORY  D9
C***TYPE      DOUBLE PRECISION (HFTI-S, DHFTI-D)
C***KEYWORDS  CURVE FITTING, LEAST SQUARES
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
C
C     This subroutine solves a linear least squares problem or a set of
C     linear least squares problems having the same matrix but different
C     right-side vectors.  The problem data consists of an M by N matrix
C     A, an M by NB matrix B, and an absolute tolerance parameter TAU
C     whose usage is described below.  The NB column vectors of B
C     represent right-side vectors for NB distinct linear least squares
C     problems.
C
C     This set of problems can also be written as the matrix least
C     squares problem
C
C                       AX = B,
C
C     where X is the N by NB solution matrix.
C
C     Note that if B is the M by M identity matrix, then X will be the
C     pseudo-inverse of A.
C
C     This subroutine first transforms the augmented matrix (A B) to a
C     matrix (R C) using premultiplying Householder transformations with
C     column interchanges.  All subdiagonal elements in the matrix R are
C     zero and its diagonal elements satisfy
C
C                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)),
C
C                       I = 1,...,L-1, where
C
C                       L = MIN(M,N).
C
C     The subroutine will compute an integer, KRANK, equal to the number
C     of diagonal terms of R that exceed TAU in magnitude. Then a
C     solution of minimum Euclidean length is computed using the first
C     KRANK rows of (R C).
C
C     To be specific we suggest that the user consider an easily
C     computable matrix norm, such as, the maximum of all column sums of
C     magnitudes.
C
C     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
C     norm of B), it is suggested that TAU be set approximately equal to
C     EPS*(norm of A).
C
C     The user must dimension all arrays appearing in the call list..
C     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
C     permits the solution of a range of problems in the same array
C     space.
C
C     The entire set of parameters for DHFTI are
C
C     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
C                       matrix A of the least squares problem AX = B.
C                       The first dimensioning parameter of the array
C                       A(*,*) is MDA, which must satisfy MDA.GE.M
C                       Either M.GE.N or M.LT.N is permitted.  There
C                       is no restriction on the rank of A.  The
C                       condition MDA.LT.M is considered an error.
C
C     B(*),MDB,NB       If NB = 0 the subroutine will perform the
C                       orthogonal decomposition but will make no
C                       references to the array B(*).  If NB.GT.0
C                       the array B(*) must initially contain the M by
C                       NB matrix B of the least squares problem AX =
C                       B.  If NB.GE.2 the array B(*) must be doubly
C                       subscripted with first dimensioning parameter
C                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may
C                       be either doubly or singly subscripted.  In
C                       the latter case the value of MDB is arbitrary
C                       but it should be set to some valid integer
C                       value such as MDB = M.
C
C                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N)
C                       is considered an error.
C
C     TAU               Absolute tolerance parameter provided by user
C                       for pseudorank determination.
C
C     H(*),G(*),IP(*)   Arrays of working space used by DHFTI.
C
C     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
C
C     A(*,*)            The contents of the array A(*,*) will be
C                       modified by the subroutine. These contents
C                       are not generally required by the user.
C
C     B(*)              On return the array B(*) will contain the N by
C                       NB solution matrix X.
C
C     KRANK             Set by the subroutine to indicate the
C                       pseudorank of A.
C
C     RNORM(*)          On return, RNORM(J) will contain the Euclidean
C                       norm of the residual vector for the problem
C                       defined by the J-th column vector of the array
C                       B(*,*) for J = 1,...,NB.
C
C     H(*),G(*)         On return these arrays respectively contain
C                       elements of the pre- and post-multiplying
C                       Householder transformations used to compute
C                       the minimum Euclidean length solution.
C
C     IP(*)             Array in which the subroutine records indices
C                       describing the permutation of column vectors.
C                       The contents of arrays H(*),G(*) and IP(*)
C                       are not generally required by the user.
C
C***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
C                 Problems, Prentice-Hall, Inc., 1974, Chapter 14.
C***ROUTINES CALLED  D1MACH, DH12, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DHFTI
      INTEGER I, II, IOPT, IP(*), IP1, J, JB, JJ, K, KP1, KRANK, L,
     *     LDIAG, LMAX, M, MDA, MDB, N, NB, NERR
      DOUBLE PRECISION A, B, D1MACH, DZERO, FACTOR,
     *     G, H, HMAX, RELEPS, RNORM, SM, SM1, SZERO, TAU, TMP
      DIMENSION A(MDA,*),B(MDB,*),H(*),G(*),RNORM(*)
      SAVE RELEPS
      DATA RELEPS /0.D0/
C     BEGIN BLOCK PERMITTING ...EXITS TO 360
C***FIRST EXECUTABLE STATEMENT  DHFTI
         IF (RELEPS.EQ.0.D0) RELEPS = D1MACH(4)
         SZERO = 0.0D0
         DZERO = 0.0D0
         FACTOR = 0.001D0
C
         K = 0
         LDIAG = MIN(M,N)
         IF (LDIAG .LE. 0) GO TO 350
C           BEGIN BLOCK PERMITTING ...EXITS TO 130
C              BEGIN BLOCK PERMITTING ...EXITS TO 120
                  IF (MDA .GE. M) GO TO 10
                     NERR = 1
                     IOPT = 2
                     CALL XERMSG ('SLATEC', 'DHFTI',
     +                  'MDA.LT.M, PROBABLE ERROR.',
     +                  NERR, IOPT)
C     ...............EXIT
                     GO TO 360
   10             CONTINUE
C
                  IF (NB .LE. 1 .OR. MAX(M,N) .LE. MDB) GO TO 20
                     NERR = 2
                     IOPT = 2
                     CALL XERMSG ('SLATEC', 'DHFTI',
     +                  'MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.',
     +                  NERR, IOPT)
C     ...............EXIT
                     GO TO 360
   20             CONTINUE
C
                  DO 100 J = 1, LDIAG
C                    BEGIN BLOCK PERMITTING ...EXITS TO 70
                        IF (J .EQ. 1) GO TO 40
C
C                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C                          ..
                           LMAX = J
                           DO 30 L = J, N
                              H(L) = H(L) - A(J-1,L)**2
                              IF (H(L) .GT. H(LMAX)) LMAX = L
   30                      CONTINUE
C                    ......EXIT
                           IF (FACTOR*H(LMAX) .GT. HMAX*RELEPS) GO TO 70
   40                   CONTINUE
C
C                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C                       ..
                        LMAX = J
                        DO 60 L = J, N
                           H(L) = 0.0D0
                           DO 50 I = J, M
                              H(L) = H(L) + A(I,L)**2
   50                      CONTINUE
                           IF (H(L) .GT. H(LMAX)) LMAX = L
   60                   CONTINUE
                        HMAX = H(LMAX)
   70                CONTINUE
C                    ..
C                     LMAX HAS BEEN DETERMINED
C
C                     DO COLUMN INTERCHANGES IF NEEDED.
C                    ..
                     IP(J) = LMAX
                     IF (IP(J) .EQ. J) GO TO 90
                        DO 80 I = 1, M
                           TMP = A(I,J)
                           A(I,J) = A(I,LMAX)
                           A(I,LMAX) = TMP
   80                   CONTINUE
                        H(LMAX) = H(J)
   90                CONTINUE
C
C                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
C                     AND B.
C                    ..
                     CALL DH12(1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,
     *                         N-J)
                     CALL DH12(2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
  100             CONTINUE
C
C                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE,
C                  TAU.
C                 ..
                  DO 110 J = 1, LDIAG
C              ......EXIT
                     IF (ABS(A(J,J)) .LE. TAU) GO TO 120
  110             CONTINUE
                  K = LDIAG
C           ......EXIT
                  GO TO 130
  120          CONTINUE
               K = J - 1
  130       CONTINUE
            KP1 = K + 1
C
C           COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
            IF (NB .LT. 1) GO TO 170
            DO 160 JB = 1, NB
               TMP = SZERO
               IF (M .LT. KP1) GO TO 150
               DO 140 I = KP1, M
                  TMP = TMP + B(I,JB)**2
  140          CONTINUE
  150          CONTINUE
               RNORM(JB) = SQRT(TMP)
  160       CONTINUE
  170       CONTINUE
C           SPECIAL FOR PSEUDORANK = 0
            IF (K .GT. 0) GO TO 210
               IF (NB .LT. 1) GO TO 200
               DO 190 JB = 1, NB
                  DO 180 I = 1, N
                     B(I,JB) = SZERO
  180             CONTINUE
  190          CONTINUE
  200          CONTINUE
            GO TO 340
  210       CONTINUE
C
C               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
C               DECOMPOSITION OF FIRST K ROWS.
C              ..
               IF (K .EQ. N) GO TO 230
                  DO 220 II = 1, K
                     I = KP1 - II
                     CALL DH12(1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  220             CONTINUE
  230          CONTINUE
C
C
               IF (NB .LT. 1) GO TO 330
               DO 320 JB = 1, NB
C
C                  SOLVE THE K BY K TRIANGULAR SYSTEM.
C                 ..
                  DO 260 L = 1, K
                     SM = DZERO
                     I = KP1 - L
                     IP1 = I + 1
                     IF (K .LT. IP1) GO TO 250
                     DO 240 J = IP1, K
                        SM = SM + A(I,J)*B(J,JB)
  240                CONTINUE
  250                CONTINUE
                     SM1 = SM
                     B(I,JB) = (B(I,JB) - SM1)/A(I,I)
  260             CONTINUE
C
C                  COMPLETE COMPUTATION OF SOLUTION VECTOR.
C                 ..
                  IF (K .EQ. N) GO TO 290
                     DO 270 J = KP1, N
                        B(J,JB) = SZERO
  270                CONTINUE
                     DO 280 I = 1, K
                        CALL DH12(2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,
     *                            MDB,1)
  280                CONTINUE
  290             CONTINUE
C
C                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C                   COLUMN INTERCHANGES.
C                 ..
                  DO 310 JJ = 1, LDIAG
                     J = LDIAG + 1 - JJ
                     IF (IP(J) .EQ. J) GO TO 300
                        L = IP(J)
                        TMP = B(L,JB)
                        B(L,JB) = B(J,JB)
                        B(J,JB) = TMP
  300                CONTINUE
  310             CONTINUE
  320          CONTINUE
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
C        ..
C         THE SOLUTION VECTORS, X, ARE NOW
C         IN THE FIRST  N  ROWS OF THE ARRAY B(,).
C
         KRANK = K
  360 CONTINUE
      RETURN
      END
