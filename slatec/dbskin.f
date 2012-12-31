*DECK DBSKIN
      SUBROUTINE DBSKIN (X, N, KODE, M, Y, NZ, IERR)
C***BEGIN PROLOGUE  DBSKIN
C***PURPOSE  Compute repeated integrals of the K-zero Bessel function.
C***LIBRARY   SLATEC
C***CATEGORY  C10F
C***TYPE      DOUBLE PRECISION (BSKIN-S, DBSKIN-D)
C***KEYWORDS  BICKLEY FUNCTIONS, EXPONENTIAL INTEGRAL,
C             INTEGRALS OF BESSEL FUNCTIONS, K-ZERO BESSEL FUNCTION
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C         The following definitions are used in DBSKIN:
C
C     Definition 1
C         KI(0,X) = K-zero Bessel function.
C
C     Definition 2
C         KI(N,X) = Bickley Function
C                 =  integral from X to infinity of KI(N-1,t)dt
C                     for X .ge. 0 and N = 1,2,...
C  _____________________________________________________________________
C    DBSKIN computes a sequence of Bickley functions (repeated integrals
C    of the K0 Bessel function); i.e. for fixed X and N and for K=1,...,
C    DBSKIN computes the sequence
C
C                     Y(K) =         KI(N+K-1,X) for KODE=1
C          or
C                     Y(K) = EXP(X)*KI(N+K-1,X) for KODE=2,
C
C         for N.ge.0 and X.ge.0 (N and X cannot be zero simultaneously).
C
C      INPUT      X is DOUBLE PRECISION
C        X      - Argument, X .ge. 0.0D0
C        N      - Order of first member of the sequence N .ge. 0
C        KODE   - Selection parameter
C             KODE = 1 returns Y(K)=        KI(N+K-1,X), K=1,M
C                  = 2 returns Y(K)=EXP(X)*KI(N+K-1,X), K=1,M
C        M      - Number of members in the sequence, M.ge.1
C
C       OUTPUT     Y is a DOUBLE PRECISION VECTOR
C         Y      - A vector of dimension at least M containing the
C                  sequence selected by KODE.
C         NZ     - Underflow flag
C                  NZ = 0 means computation completed
C                     = 1 means an exponential underflow occurred on
C                         KODE=1.  Y(K)=0.0D0, K=1,...,M is returned
C                         KODE=1 AND Y(K)=0.0E0, K=1,...,M IS RETURNED
C         IERR   - Error flag
C                    IERR=0, Normal return, computation completed
C                    IERR=1, Input error,   no computation
C                    IERR=2, Error,         no computation
C                            Algorithm termination condition not met
C
C         The nominal computational accuracy is the maximum of unit
C         roundoff (=D1MACH(4)) and 1.0D-18 since critical constants
C         are given to only 18 digits.
C
C         BSKIN is the single precision version of DBSKIN.
C
C *Long Description:
C
C         Numerical recurrence on
C
C      (L-1)*KI(L,X) = X(KI(L-3,X) - KI(L-1,X)) + (L-2)*KI(L-2,X)
C
C         is stable where recurrence is carried forward or backward
C         away from INT(X+0.5).  The power series for indices 0,1 and 2
C         on 0.le.X.le.2 starts a stable recurrence for indices
C         greater than 2.  If N is sufficiently large (N.gt.NLIM), the
C         uniform asymptotic expansion for N to INFINITY is more
C         economical.  On X.gt.2 the recursion is started by evaluating
C         the uniform expansion for the three members whose indices are
C         closest to INT(X+0.5) within the set N,...,N+M-1.  Forward
C         recurrence, backward recurrence or both complete the
C         sequence depending on the relation of INT(X+0.5) to the
C         indices N,...,N+M-1.
C
C***REFERENCES  D. E. Amos, Uniform asymptotic expansions for
C                 exponential integrals E(N,X) and Bickley functions
C                 KI(N,X), ACM Transactions on Mathematical Software,
C                 1983.
C               D. E. Amos, A portable Fortran subroutine for the
C                 Bickley functions KI(N,X), Algorithm 609, ACM
C                 Transactions on Mathematical Software, 1983.
C***ROUTINES CALLED  D1MACH, DBKIAS, DBKISR, DEXINT, DGAMRN, I1MACH
C***REVISION HISTORY  (YYMMDD)
C   820601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891009  Removed unreferenced statement label.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBSKIN
      INTEGER I, ICASE, IERR, IL, I1M, K, KK, KODE, KTRMS, M,
     * M3, N, NE, NFLG, NL, NLIM, NN, NP, NS, NT, NZ
      INTEGER I1MACH
      DOUBLE PRECISION A, ENLIM, EXI, FN, GR, H, HN, HRTPI, SS, TOL,
     * T1, T2, W, X, XLIM, XNLIM, XP, Y, YS, YSS
      DOUBLE PRECISION DGAMRN, D1MACH
      DIMENSION EXI(102), A(50), YS(3), YSS(3), H(31), Y(*)
      SAVE A, HRTPI
C-----------------------------------------------------------------------
C             COEFFICIENTS IN SERIES OF EXPONENTIAL INTEGRALS
C-----------------------------------------------------------------------
      DATA A(1), A(2), A(3), A(4), A(5), A(6), A(7), A(8), A(9), A(10),
     * A(11), A(12), A(13), A(14), A(15), A(16), A(17), A(18), A(19),
     * A(20), A(21), A(22), A(23), A(24) /1.00000000000000000D+00,
     * 5.00000000000000000D-01,3.75000000000000000D-01,
     * 3.12500000000000000D-01,2.73437500000000000D-01,
     * 2.46093750000000000D-01,2.25585937500000000D-01,
     * 2.09472656250000000D-01,1.96380615234375000D-01,
     * 1.85470581054687500D-01,1.76197052001953125D-01,
     * 1.68188095092773438D-01,1.61180257797241211D-01,
     * 1.54981017112731934D-01,1.49445980787277222D-01,
     * 1.44464448094367981D-01,1.39949934091418982D-01,
     * 1.35833759559318423D-01,1.32060599571559578D-01,
     * 1.28585320635465905D-01,1.25370687619579257D-01,
     * 1.22385671247684513D-01,1.19604178719328047D-01,
     * 1.17004087877603524D-01/
      DATA A(25), A(26), A(27), A(28), A(29), A(30), A(31), A(32),
     * A(33), A(34), A(35), A(36), A(37), A(38), A(39), A(40), A(41),
     * A(42), A(43), A(44), A(45), A(46), A(47), A(48)
     * /1.14566502713486784D-01,1.12275172659217048D-01,
     * 1.10116034723462874D-01,1.08076848895250599D-01,
     * 1.06146905164978267D-01,1.04316786110409676D-01,
     * 1.02578173008569515D-01,1.00923686347140974D-01,
     * 9.93467537479668965D-02,9.78414999033007314D-02,
     * 9.64026543164874854D-02,9.50254735405376642D-02,
     * 9.37056752969190855D-02,9.24393823875012600D-02,
     * 9.12230747245078224D-02,9.00535481254756708D-02,
     * 8.89278787739072249D-02,8.78433924473961612D-02,
     * 8.67976377754033498D-02,8.57883629175498224D-02,
     * 8.48134951571231199D-02,8.38711229887106408D-02,
     * 8.29594803475290034D-02,8.20769326842574183D-02/
      DATA A(49), A(50) /8.12219646354630702D-02,8.03931690779583449D-02
     * /
C-----------------------------------------------------------------------
C             SQRT(PI)/2
C-----------------------------------------------------------------------
      DATA HRTPI /8.86226925452758014D-01/
C
C***FIRST EXECUTABLE STATEMENT  DBSKIN
      IERR = 0
      NZ=0
      IF (X.LT.0.0D0) IERR=1
      IF (N.LT.0) IERR=1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR=1
      IF (M.LT.1) IERR=1
      IF (X.EQ.0.0D0 .AND. N.EQ.0) IERR=1
      IF (IERR.NE.0) RETURN
      IF (X.EQ.0.0D0) GO TO 300
      I1M = -I1MACH(15)
      T1 = 2.3026D0*D1MACH(5)*I1M
      XLIM = T1 - 3.228086D0
      T2 = T1 + (N+M-1)
      IF (T2.GT.1000.0D0) XLIM = T1 - 0.5D0*(LOG(T2)-0.451583D0)
      IF (X.GT.XLIM .AND. KODE.EQ.1) GO TO 320
      TOL = MAX(D1MACH(4),1.0D-18)
      I1M = I1MACH(14)
C-----------------------------------------------------------------------
C     LN(NLIM) = 0.125*LN(EPS),   NLIM = 2*KTRMS+N
C-----------------------------------------------------------------------
      XNLIM = 0.287823D0*(I1M-1)*D1MACH(5)
      ENLIM = EXP(XNLIM)
      NLIM = INT(ENLIM) + 2
      NLIM = MIN(100,NLIM)
      NLIM = MAX(20,NLIM)
      M3 = MIN(M,3)
      NL = N + M - 1
      IF (X.GT.2.0D0) GO TO 130
      IF (N.GT.NLIM) GO TO 280
C-----------------------------------------------------------------------
C     COMPUTATION BY SERIES FOR 0.LE.X.LE.2
C-----------------------------------------------------------------------
      NFLG = 0
      NN = N
      IF (NL.LE.2) GO TO 60
      M3 = 3
      NN = 0
      NFLG = 1
   60 CONTINUE
      XP = 1.0D0
      IF (KODE.EQ.2) XP = EXP(X)
      DO 80 I=1,M3
        CALL DBKISR(X, NN, W, IERR)
      IF(IERR.NE.0) RETURN
        W = W*XP
        IF (NN.LT.N) GO TO 70
        KK = NN - N + 1
        Y(KK) = W
   70   CONTINUE
        YS(I) = W
        NN = NN + 1
   80 CONTINUE
      IF (NFLG.EQ.0) RETURN
      NS = NN
      XP = 1.0D0
   90 CONTINUE
C-----------------------------------------------------------------------
C     FORWARD RECURSION SCALED BY EXP(X) ON ICASE=0,1,2
C-----------------------------------------------------------------------
      FN = NS - 1
      IL = NL - NS + 1
      IF (IL.LE.0) RETURN
      DO 110 I=1,IL
        T1 = YS(2)
        T2 = YS(3)
        YS(3) = (X*(YS(1)-YS(3))+(FN-1.0D0)*YS(2))/FN
        YS(2) = T2
        YS(1) = T1
        FN = FN + 1.0D0
        IF (NS.LT.N) GO TO 100
        KK = NS - N + 1
        Y(KK) = YS(3)*XP
  100   CONTINUE
        NS = NS + 1
  110 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     COMPUTATION BY ASYMPTOTIC EXPANSION FOR X.GT.2
C-----------------------------------------------------------------------
  130 CONTINUE
      W = X + 0.5D0
      NT = INT(W)
      IF (NL.GT.NT) GO TO 270
C-----------------------------------------------------------------------
C     CASE NL.LE.NT, ICASE=0
C-----------------------------------------------------------------------
      ICASE = 0
      NN = NL
      NFLG = MIN(M-M3,1)
  140 CONTINUE
      KK = (NLIM-NN)/2
      KTRMS = MAX(0,KK)
      NS = NN + 1
      NP = NN - M3 + 1
      XP = 1.0D0
      IF (KODE.EQ.1) XP = EXP(-X)
      DO 150 I=1,M3
        KK = I
        CALL DBKIAS(X, NP, KTRMS, A, W, KK, NE, GR, H, IERR)
      IF(IERR.NE.0) RETURN
        YS(I) = W
        NP = NP + 1
  150 CONTINUE
C-----------------------------------------------------------------------
C     SUM SERIES OF EXPONENTIAL INTEGRALS BACKWARD
C-----------------------------------------------------------------------
      IF (KTRMS.EQ.0) GO TO 160
      NE = KTRMS + KTRMS + 1
      NP = NN - M3 + 2
      CALL DEXINT(X, NP, 2, NE, TOL, EXI, NZ, IERR)
      IF (NZ.NE.0) GO TO 320
  160 CONTINUE
      DO 190 I=1,M3
        SS = 0.0D0
        IF (KTRMS.EQ.0) GO TO 180
        KK = I + KTRMS + KTRMS - 2
        IL = KTRMS
        DO 170 K=1,KTRMS
          SS = SS + A(IL)*EXI(KK)
          KK = KK - 2
          IL = IL - 1
  170   CONTINUE
  180   CONTINUE
        YS(I) = YS(I) + SS
  190 CONTINUE
      IF (ICASE.EQ.1) GO TO 200
      IF (NFLG.NE.0) GO TO 220
  200 CONTINUE
      DO 210 I=1,M3
        Y(I) = YS(I)*XP
  210 CONTINUE
      IF (ICASE.EQ.1 .AND. NFLG.EQ.1) GO TO 90
      RETURN
  220 CONTINUE
C-----------------------------------------------------------------------
C     BACKWARD RECURSION SCALED BY EXP(X) ICASE=0,2
C-----------------------------------------------------------------------
      KK = NN - N + 1
      K = M3
      DO 230 I=1,M3
        Y(KK) = YS(K)*XP
        YSS(I) = YS(I)
        KK = KK - 1
        K = K - 1
  230 CONTINUE
      IL = KK
      IF (IL.LE.0) GO TO 250
      FN = NN - 3
      DO 240 I=1,IL
        T1 = YS(2)
        T2 = YS(1)
        YS(1) = YS(2) + ((FN+2.0D0)*YS(3)-(FN+1.0D0)*YS(1))/X
        YS(2) = T2
        YS(3) = T1
        Y(KK) = YS(1)*XP
        KK = KK - 1
        FN = FN - 1.0D0
  240 CONTINUE
  250 CONTINUE
      IF (ICASE.NE.2) RETURN
      DO 260 I=1,M3
        YS(I) = YSS(I)
  260 CONTINUE
      GO TO 90
  270 CONTINUE
      IF (N.LT.NT) GO TO 290
C-----------------------------------------------------------------------
C     ICASE=1, NT.LE.N.LE.NL WITH FORWARD RECURSION
C-----------------------------------------------------------------------
  280 CONTINUE
      NN = N + M3 - 1
      NFLG = MIN(M-M3,1)
      ICASE = 1
      GO TO 140
C-----------------------------------------------------------------------
C     ICASE=2, N.LT.NT.LT.NL WITH BOTH FORWARD AND BACKWARD RECURSION
C-----------------------------------------------------------------------
  290 CONTINUE
      NN = NT + 1
      NFLG = MIN(M-M3,1)
      ICASE = 2
      GO TO 140
C-----------------------------------------------------------------------
C     X=0 CASE
C-----------------------------------------------------------------------
  300 CONTINUE
      FN = N
      HN = 0.5D0*FN
      GR = DGAMRN(HN)
      Y(1) = HRTPI*GR
      IF (M.EQ.1) RETURN
      Y(2) = HRTPI/(HN*GR)
      IF (M.EQ.2) RETURN
      DO 310 K=3,M
        Y(K) = FN*Y(K-2)/(FN+1.0D0)
        FN = FN + 1.0D0
  310 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     UNDERFLOW ON KODE=1, X.GT.XLIM
C-----------------------------------------------------------------------
  320 CONTINUE
      NZ=M
      DO 330 I=1,M
        Y(I) = 0.0D0
  330 CONTINUE
      RETURN
      END
