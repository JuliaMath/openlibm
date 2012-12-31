*DECK DMGSBV
      SUBROUTINE DMGSBV (M, N, A, IA, NIV, IFLAG, S, P, IP, INHOMO, V,
     +   W, WCND)
C***BEGIN PROLOGUE  DMGSBV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBVSUP
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (MGSBV-S, DMGSBV-D)
C***AUTHOR  Watts, H. A., (SNLA)
C***DESCRIPTION
C
C **********************************************************************
C Orthogonalize a set of N double precision vectors and determine their
C rank.
C
C **********************************************************************
C INPUT
C **********************************************************************
C   M = dimension of vectors.
C   N = no. of vectors.
C   A = array whose first N cols contain the vectors.
C   IA = first dimension of array A (col length).
C   NIV = number of independent vectors needed.
C   INHOMO = 1 corresponds to having a non-zero particular solution.
C   V = particular solution vector (not included in the pivoting).
C   INDPVT = 1 means pivoting will not be used.
C
C **********************************************************************
C OUTPUT
C **********************************************************************
C   NIV = no. of linear independent vectors in input set.
C     A = matrix whose first NIV cols. contain NIV orthogonal vectors
C         which span the vector space determined by the input vectors.
C   IFLAG
C          = 0 success
C          = 1 incorrect input
C          = 2 rank of new vectors less than N
C   P = decomposition matrix.  P is upper triangular and
C             (old vectors) = (new vectors) * P.
C         The old vectors will be reordered due to pivoting.
C         The dimension of P must be .GE. N*(N+1)/2.
C             (  N*(2*N+1) when N .NE. NFCC )
C   IP = pivoting vector. The dimension of IP must be .GE. N.
C             (  2*N when N .NE. NFCC )
C   S = square of norms of incoming vectors.
C   V = vector which is orthogonal to the vectors of A.
C   W = orthogonalization information for the vector V.
C   WCND = worst case (smallest) norm decrement value of the
C          vectors being orthogonalized  (represents a test
C          for linear dependence of the vectors).
C **********************************************************************
C
C***SEE ALSO  DBVSUP
C***ROUTINES CALLED  DDOT, DPRVEC
C***COMMON BLOCKS    DML18J, DML5MC
C***REVISION HISTORY  (YYMMDD)
C   750601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890921  Realigned order of variables in certain COMMON blocks.
C           (WRB)
C   890921  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C   910722  Updated AUTHOR section.  (ALS)
C***END PROLOGUE  DMGSBV
C
      DOUBLE PRECISION DDOT, DPRVEC
      INTEGER I, IA, ICOCO, IFLAG, INDPVT, INHOMO, INTEG, IP(*), IP1,
     1     IX, IZ, J, JK, JP, JQ, JY, JZ, K, KD, KJ, KP, L, LIX, LPAR,
     2     LR, M, M2, MXNON, N, NDISK, NEQ, NEQIVP, NFCC, NIC, NIV,
     3     NIVN, NMNR, NN, NOPG, NP1, NPS, NR, NRM1, NTAPE, NTP,
     4     NUMORT, NXPTS
      DOUBLE PRECISION A(IA,*), AE, DOT, EPS, FOURU, P(*), PJP, PSAVE,
     1     RE, RY, S(*), SQOVFL, SRU, SV, T, TOL, TWOU, URO, V(*), VL,
     2     VNORM, W(*), WCND, Y
C
C
      COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ,
     1                INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC,
     2                ICOCO
C
      COMMON /DML5MC/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
C
C***FIRST EXECUTABLE STATEMENT  DMGSBV
      IF (M .GT. 0 .AND. N .GT. 0 .AND. IA .GE. M) GO TO 10
         IFLAG = 1
      GO TO 280
   10 CONTINUE
C        BEGIN BLOCK PERMITTING ...EXITS TO 270
C           BEGIN BLOCK PERMITTING ...EXITS TO 260
C
               JP = 0
               IFLAG = 0
               NP1 = N + 1
               Y = 0.0D0
               M2 = M/2
C
C              CALCULATE SQUARE OF NORMS OF INCOMING VECTORS AND SEARCH
C              FOR VECTOR WITH LARGEST MAGNITUDE
C
               J = 0
               DO 40 I = 1, N
                  VL = DDOT(M,A(1,I),1,A(1,I),1)
                  S(I) = VL
                  IF (N .EQ. NFCC) GO TO 20
                     J = 2*I - 1
                     P(J) = VL
                     IP(J) = J
   20             CONTINUE
                  J = J + 1
                  P(J) = VL
                  IP(J) = J
                  IF (VL .LE. Y) GO TO 30
                     Y = VL
                     IX = I
   30             CONTINUE
   40          CONTINUE
               IF (INDPVT .NE. 1) GO TO 50
                  IX = 1
                  Y = P(1)
   50          CONTINUE
               LIX = IX
               IF (N .NE. NFCC) LIX = 2*IX - 1
               P(LIX) = P(1)
               S(NP1) = 0.0D0
               IF (INHOMO .EQ. 1) S(NP1) = DDOT(M,V,1,V,1)
               WCND = 1.0D0
               NIVN = NIV
               NIV = 0
C
C           ...EXIT
               IF (Y .EQ. 0.0D0) GO TO 260
C              *********************************************************
               DO 240 NR = 1, N
C                 BEGIN BLOCK PERMITTING ...EXITS TO 230
C              ......EXIT
                     IF (NIVN .EQ. NIV) GO TO 250
                     NIV = NR
                     IF (IX .EQ. NR) GO TO 130
C
C                       PIVOTING OF COLUMNS OF P MATRIX
C
                        NN = N
                        LIX = IX
                        LR = NR
                        IF (N .EQ. NFCC) GO TO 60
                           NN = NFCC
                           LIX = 2*IX - 1
                           LR = 2*NR - 1
   60                   CONTINUE
                        IF (NR .EQ. 1) GO TO 80
                           KD = LIX - LR
                           KJ = LR
                           NRM1 = LR - 1
                           DO 70 J = 1, NRM1
                              PSAVE = P(KJ)
                              JK = KJ + KD
                              P(KJ) = P(JK)
                              P(JK) = PSAVE
                              KJ = KJ + NN - J
   70                      CONTINUE
                           JY = JK + NMNR
                           JZ = JY - KD
                           P(JY) = P(JZ)
   80                   CONTINUE
                        IZ = IP(LIX)
                        IP(LIX) = IP(LR)
                        IP(LR) = IZ
                        SV = S(IX)
                        S(IX) = S(NR)
                        S(NR) = SV
                        IF (N .EQ. NFCC) GO TO 110
                           IF (NR .EQ. 1) GO TO 100
                              KJ = LR + 1
                              DO 90 K = 1, NRM1
                                 PSAVE = P(KJ)
                                 JK = KJ + KD
                                 P(KJ) = P(JK)
                                 P(JK) = PSAVE
                                 KJ = KJ + NFCC - K
   90                         CONTINUE
  100                      CONTINUE
                           IZ = IP(LIX+1)
                           IP(LIX+1) = IP(LR+1)
                           IP(LR+1) = IZ
  110                   CONTINUE
C
C                       PIVOTING OF COLUMNS OF VECTORS
C
                        DO 120 L = 1, M
                           T = A(L,IX)
                           A(L,IX) = A(L,NR)
                           A(L,NR) = T
  120                   CONTINUE
  130                CONTINUE
C
C                    CALCULATE P(NR,NR) AS NORM SQUARED OF PIVOTAL
C                    VECTOR
C
                     JP = JP + 1
                     P(JP) = Y
                     RY = 1.0D0/Y
                     NMNR = N - NR
                     IF (N .EQ. NFCC) GO TO 140
                        NMNR = NFCC - (2*NR - 1)
                        JP = JP + 1
                        P(JP) = 0.0D0
                        KP = JP + NMNR
                        P(KP) = Y
  140                CONTINUE
                     IF (NR .EQ. N .OR. NIVN .EQ. NIV) GO TO 200
C
C                       CALCULATE ORTHOGONAL PROJECTION VECTORS AND
C                       SEARCH FOR LARGEST NORM
C
                        Y = 0.0D0
                        IP1 = NR + 1
                        IX = IP1
C                       ************************************************
                        DO 190 J = IP1, N
                           DOT = DDOT(M,A(1,NR),1,A(1,J),1)
                           JP = JP + 1
                           JQ = JP + NMNR
                           IF (N .NE. NFCC) JQ = JQ + NMNR - 1
                           P(JQ) = P(JP) - DOT*(DOT*RY)
                           P(JP) = DOT*RY
                           DO 150 I = 1, M
                              A(I,J) = A(I,J) - P(JP)*A(I,NR)
  150                      CONTINUE
                           IF (N .EQ. NFCC) GO TO 170
                              KP = JP + NMNR
                              JP = JP + 1
                              PJP = RY*DPRVEC(M,A(1,NR),A(1,J))
                              P(JP) = PJP
                              P(KP) = -PJP
                              KP = KP + 1
                              P(KP) = RY*DOT
                              DO 160 K = 1, M2
                                 L = M2 + K
                                 A(K,J) = A(K,J) - PJP*A(L,NR)
                                 A(L,J) = A(L,J) + PJP*A(K,NR)
  160                         CONTINUE
                              P(JQ) = P(JQ) - PJP*(PJP/RY)
  170                      CONTINUE
C
C                          TEST FOR CANCELLATION IN RECURRENCE RELATION
C
                           IF (P(JQ) .LE. S(J)*SRU)
     1                        P(JQ) = DDOT(M,A(1,J),1,A(1,J),1)
                           IF (P(JQ) .LE. Y) GO TO 180
                              Y = P(JQ)
                              IX = J
  180                      CONTINUE
  190                   CONTINUE
                        IF (N .NE. NFCC) JP = KP
C                       ************************************************
                        IF (INDPVT .EQ. 1) IX = IP1
C
C                       RECOMPUTE NORM SQUARED OF PIVOTAL VECTOR WITH
C                       SCALAR PRODUCT
C
                        Y = DDOT(M,A(1,IX),1,A(1,IX),1)
C           ............EXIT
                        IF (Y .LE. EPS*S(IX)) GO TO 260
                        WCND = MIN(WCND,Y/S(IX))
  200                CONTINUE
C
C                    COMPUTE ORTHOGONAL PROJECTION OF PARTICULAR
C                    SOLUTION
C
C                 ...EXIT
                     IF (INHOMO .NE. 1) GO TO 230
                     LR = NR
                     IF (N .NE. NFCC) LR = 2*NR - 1
                     W(LR) = DDOT(M,A(1,NR),1,V,1)*RY
                     DO 210 I = 1, M
                        V(I) = V(I) - W(LR)*A(I,NR)
  210                CONTINUE
C                 ...EXIT
                     IF (N .EQ. NFCC) GO TO 230
                     LR = 2*NR
                     W(LR) = RY*DPRVEC(M,V,A(1,NR))
                     DO 220 K = 1, M2
                        L = M2 + K
                        V(K) = V(K) + W(LR)*A(L,NR)
                        V(L) = V(L) - W(LR)*A(K,NR)
  220                CONTINUE
  230             CONTINUE
  240          CONTINUE
  250          CONTINUE
C              *********************************************************
C
C                  TEST FOR LINEAR DEPENDENCE OF PARTICULAR SOLUTION
C
C        ......EXIT
               IF (INHOMO .NE. 1) GO TO 270
               IF ((N .GT. 1) .AND. (S(NP1) .LT. 1.0)) GO TO 270
               VNORM = DDOT(M,V,1,V,1)
               IF (S(NP1) .NE. 0.0D0) WCND = MIN(WCND,VNORM/S(NP1))
C        ......EXIT
               IF (VNORM .GE. EPS*S(NP1)) GO TO 270
  260       CONTINUE
            IFLAG = 2
            WCND = EPS
  270    CONTINUE
  280 CONTINUE
      RETURN
      END
