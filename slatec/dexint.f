*DECK DEXINT
      SUBROUTINE DEXINT (X, N, KODE, M, TOL, EN, NZ, IERR)
C***BEGIN PROLOGUE  DEXINT
C***PURPOSE  Compute an M member sequence of exponential integrals
C            E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.
C***LIBRARY   SLATEC
C***CATEGORY  C5
C***TYPE      DOUBLE PRECISION (EXINT-S, DEXINT-D)
C***KEYWORDS  EXPONENTIAL INTEGRAL, SPECIAL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C         DEXINT computes M member sequences of exponential integrals
C         E(N+K,X), K=0,1,...,M-1 for N .GE. 1 and X .GE. 0.  The
C         exponential integral is defined by
C
C         E(N,X)=integral on (1,infinity) of EXP(-XT)/T**N
C
C         where X=0.0 and N=1 cannot occur simultaneously.  Formulas
C         and notation are found in the NBS Handbook of Mathematical
C         Functions (ref. 1).
C
C         The power series is implemented for X .LE. XCUT and the
C         confluent hypergeometric representation
C
C                     E(A,X) = EXP(-X)*(X**(A-1))*U(A,A,X)
C
C         is computed for X .GT. XCUT.  Since sequences are computed in
C         a stable fashion by recurring away from X, A is selected as
C         the integer closest to X within the constraint N .LE. A .LE.
C         N+M-1.  For the U computation, A is further modified to be the
C         nearest even integer.  Indices are carried forward or
C         backward by the two term recursion relation
C
C                     K*E(K+1,X) + X*E(K,X) = EXP(-X)
C
C         once E(A,X) is computed.  The U function is computed by means
C         of the backward recursive Miller algorithm applied to the
C         three term contiguous relation for U(A+K,A,X), K=0,1,...
C         This produces accurate ratios and determines U(A+K,A,X), and
C         hence E(A,X), to within a multiplicative constant C.
C         Another contiguous relation applied to C*U(A,A,X) and
C         C*U(A+1,A,X) gets C*U(A+1,A+1,X), a quantity proportional to
C         E(A+1,X).  The normalizing constant C is obtained from the
C         two term recursion relation above with K=A.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C     Description of Arguments
C
C         Input     * X and TOL are double precision *
C           X       X .GT. 0.0 for N=1 and  X .GE. 0.0 for N .GE. 2
C           N       order of the first member of the sequence, N .GE. 1
C                   (X=0.0 and N=1 is an error)
C           KODE    a selection parameter for scaled values
C                   KODE=1   returns        E(N+K,X), K=0,1,...,M-1.
C                       =2   returns EXP(X)*E(N+K,X), K=0,1,...,M-1.
C           M       number of exponential integrals in the sequence,
C                   M .GE. 1
C           TOL     relative accuracy wanted, ETOL .LE. TOL .LE. 0.1
C                   ETOL is the larger of double precision unit
C                   roundoff = D1MACH(4) and 1.0D-18
C
C         Output    * EN is a double precision vector *
C           EN      a vector of dimension at least M containing values
C                   EN(K) = E(N+K-1,X) or EXP(X)*E(N+K-1,X), K=1,M
C                   depending on KODE
C           NZ      underflow indicator
C                   NZ=0   a normal return
C                   NZ=M   X exceeds XLIM and an underflow occurs.
C                          EN(K)=0.0D0 , K=1,M returned on KODE=1
C           IERR    error flag
C                   IERR=0, normal return, computation completed
C                   IERR=1, input error,   no computation
C                   IERR=2, error,         no computation
C                           algorithm termination condition not met
C
C***REFERENCES  M. Abramowitz and I. A. Stegun, Handbook of
C                 Mathematical Functions, NBS AMS Series 55, U.S. Dept.
C                 of Commerce, 1955.
C               D. E. Amos, Computation of exponential integrals, ACM
C                 Transactions on Mathematical Software 6, (1980),
C                 pp. 365-377 and pp. 420-428.
C***ROUTINES CALLED  D1MACH, DPSIXN, I1MACH
C***REVISION HISTORY  (YYMMDD)
C   800501  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   910408  Updated the REFERENCES section.  (WRB)
C   920207  Updated with code with a revision date of 880811 from
C           D. Amos.  Included correction of argument list.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DEXINT
      DOUBLE PRECISION A,AA,AAMS,AH,AK,AT,B,BK,BT,CC,CNORM,CT,EM,EMX,EN,
     1                 ETOL,FNM,FX,PT,P1,P2,S,TOL,TX,X,XCUT,XLIM,XTOL,Y,
     2                 YT,Y1,Y2
      DOUBLE PRECISION D1MACH,DPSIXN
      INTEGER I,IC,ICASE,ICT,IERR,IK,IND,IX,I1M,JSET,K,KK,KN,KODE,KS,M,
     1        ML,MU,N,ND,NM,NZ
      INTEGER I1MACH
      DIMENSION EN(*), A(99), B(99), Y(2)
      SAVE XCUT
      DATA XCUT / 2.0D0 /
C***FIRST EXECUTABLE STATEMENT  DEXINT
      IERR = 0
      NZ = 0
      ETOL = MAX(D1MACH(4),0.5D-18)
      IF (X.LT.0.0D0) IERR = 1
      IF (N.LT.1) IERR = 1
      IF (KODE.LT.1 .OR. KODE.GT.2) IERR = 1
      IF (M.LT.1) IERR = 1
      IF (TOL.LT.ETOL .OR. TOL.GT.0.1D0) IERR = 1
      IF (X.EQ.0.0D0 .AND. N.EQ.1) IERR = 1
      IF(IERR.NE.0) RETURN
      I1M = -I1MACH(15)
      PT = 2.3026D0*I1M*D1MACH(5)
      XLIM = PT - 6.907755D0
      BT = PT + (N+M-1)
      IF (BT.GT.1000.0D0) XLIM = PT - LOG(BT)
C
      IF (X.GT.XCUT) GO TO 100
      IF (X.EQ.0.0D0 .AND. N.GT.1) GO TO 80
C-----------------------------------------------------------------------
C     SERIES FOR E(N,X) FOR X.LE.XCUT
C-----------------------------------------------------------------------
      TX = X + 0.5D0
      IX = TX
C-----------------------------------------------------------------------
C     ICASE=1 MEANS INTEGER CLOSEST TO X IS 2 AND N=1
C     ICASE=2 MEANS INTEGER CLOSEST TO X IS 0,1, OR 2 AND N.GE.2
C-----------------------------------------------------------------------
      ICASE = 2
      IF (IX.GT.N) ICASE = 1
      NM = N - ICASE + 1
      ND = NM + 1
      IND = 3 - ICASE
      MU = M - IND
      ML = 1
      KS = ND
      FNM = NM
      S = 0.0D0
      XTOL = 3.0D0*TOL
      IF (ND.EQ.1) GO TO 10
      XTOL = 0.3333D0*TOL
      S = 1.0D0/FNM
   10 CONTINUE
      AA = 1.0D0
      AK = 1.0D0
      IC = 35
      IF (X.LT.ETOL) IC = 1
      DO 50 I=1,IC
        AA = -AA*X/AK
        IF (I.EQ.NM) GO TO 30
        S = S - AA/(AK-FNM)
        IF (ABS(AA).LE.XTOL*ABS(S)) GO TO 20
        AK = AK + 1.0D0
        GO TO 50
   20   CONTINUE
        IF (I.LT.2) GO TO 40
        IF (ND-2.GT.I .OR. I.GT.ND-1) GO TO 60
        AK = AK + 1.0D0
        GO TO 50
   30   S = S + AA*(-LOG(X)+DPSIXN(ND))
        XTOL = 3.0D0*TOL
   40   AK = AK + 1.0D0
   50 CONTINUE
      IF (IC.NE.1) GO TO 340
   60 IF (ND.EQ.1) S = S + (-LOG(X)+DPSIXN(1))
      IF (KODE.EQ.2) S = S*EXP(X)
      EN(1) = S
      EMX = 1.0D0
      IF (M.EQ.1) GO TO 70
      EN(IND) = S
      AA = KS
      IF (KODE.EQ.1) EMX = EXP(-X)
      GO TO (220, 240), ICASE
   70 IF (ICASE.EQ.2) RETURN
      IF (KODE.EQ.1) EMX = EXP(-X)
      EN(1) = (EMX-S)/X
      RETURN
   80 CONTINUE
      DO 90 I=1,M
        EN(I) = 1.0D0/(N+I-2)
   90 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     BACKWARD RECURSIVE MILLER ALGORITHM FOR
C              E(N,X)=EXP(-X)*(X**(N-1))*U(N,N,X)
C     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO X.
C     U(A,B,X) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
C-----------------------------------------------------------------------
  100 CONTINUE
      EMX = 1.0D0
      IF (KODE.EQ.2) GO TO 130
      IF (X.LE.XLIM) GO TO 120
      NZ = M
      DO 110 I=1,M
        EN(I) = 0.0D0
  110 CONTINUE
      RETURN
  120 EMX = EXP(-X)
  130 CONTINUE
      TX = X + 0.5D0
      IX = TX
      KN = N + M - 1
      IF (KN.LE.IX) GO TO 140
      IF (N.LT.IX .AND. IX.LT.KN) GO TO 170
      IF (N.GE.IX) GO TO 160
      GO TO 340
  140 ICASE = 1
      KS = KN
      ML = M - 1
      MU = -1
      IND = M
      IF (KN.GT.1) GO TO 180
  150 KS = 2
      ICASE = 3
      GO TO 180
  160 ICASE = 2
      IND = 1
      KS = N
      MU = M - 1
      IF (N.GT.1) GO TO 180
      IF (KN.EQ.1) GO TO 150
      IX = 2
  170 ICASE = 1
      KS = IX
      ML = IX - N
      IND = ML + 1
      MU = KN - IX
  180 CONTINUE
      IK = KS/2
      AH = IK
      JSET = 1 + KS - (IK+IK)
C-----------------------------------------------------------------------
C     START COMPUTATION FOR
C              EN(IND) = C*U( A , A ,X)    JSET=1
C              EN(IND) = C*U(A+1,A+1,X)    JSET=2
C     FOR AN EVEN INTEGER A.
C-----------------------------------------------------------------------
      IC = 0
      AA = AH + AH
      AAMS = AA - 1.0D0
      AAMS = AAMS*AAMS
      TX = X + X
      FX = TX + TX
      AK = AH
      XTOL = TOL
      IF (TOL.LE.1.0D-3) XTOL = 20.0D0*TOL
      CT = AAMS + FX*AH
      EM = (AH+1.0D0)/((X+AA)*XTOL*SQRT(CT))
      BK = AA
      CC = AH*AH
C-----------------------------------------------------------------------
C     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
C     RECURSION
C-----------------------------------------------------------------------
      P1 = 0.0D0
      P2 = 1.0D0
  190 CONTINUE
      IF (IC.EQ.99) GO TO 340
      IC = IC + 1
      AK = AK + 1.0D0
      AT = BK/(BK+AK+CC+IC)
      BK = BK + AK + AK
      A(IC) = AT
      BT = (AK+AK+X)/(AK+1.0D0)
      B(IC) = BT
      PT = P2
      P2 = BT*P2 - AT*P1
      P1 = PT
      CT = CT + FX
      EM = EM*AT*(1.0D0-TX/CT)
      IF (EM*(AK+1.0D0).GT.P1*P1) GO TO 190
      ICT = IC
      KK = IC + 1
      BT = TX/(CT+FX)
      Y2 = (BK/(BK+CC+KK))*(P1/P2)*(1.0D0-BT+0.375D0*BT*BT)
      Y1 = 1.0D0
C-----------------------------------------------------------------------
C     BACKWARD RECURRENCE FOR
C              Y1=             C*U( A ,A,X)
C              Y2= C*(A/(1+A/2))*U(A+1,A,X)
C-----------------------------------------------------------------------
      DO 200 K=1,ICT
        KK = KK - 1
        YT = Y1
        Y1 = (B(KK)*Y1-Y2)/A(KK)
        Y2 = YT
  200 CONTINUE
C-----------------------------------------------------------------------
C     THE CONTIGUOUS RELATION
C              X*U(B,C+1,X)=(C-B)*U(B,C,X)+U(B-1,C,X)
C     WITH  B=A+1 , C=A IS USED FOR
C              Y(2) = C * U(A+1,A+1,X)
C     X IS INCORPORATED INTO THE NORMALIZING RELATION
C-----------------------------------------------------------------------
      PT = Y2/Y1
      CNORM = 1.0E0 - PT*(AH+1.0E0)/AA
      Y(1) = 1.0E0/(CNORM*AA+X)
      Y(2) = CNORM*Y(1)
      IF (ICASE.EQ.3) GO TO 210
      EN(IND) = EMX*Y(JSET)
      IF (M.EQ.1) RETURN
      AA = KS
      GO TO (220, 240), ICASE
C-----------------------------------------------------------------------
C     RECURSION SECTION  N*E(N+1,X) + X*E(N,X)=EMX
C-----------------------------------------------------------------------
  210 EN(1) = EMX*(1.0E0-Y(1))/X
      RETURN
  220 K = IND - 1
      DO 230 I=1,ML
        AA = AA - 1.0D0
        EN(K) = (EMX-AA*EN(K+1))/X
        K = K - 1
  230 CONTINUE
      IF (MU.LE.0) RETURN
      AA = KS
  240 K = IND
      DO 250 I=1,MU
        EN(K+1) = (EMX-X*EN(K))/AA
        AA = AA + 1.0D0
        K = K + 1
  250 CONTINUE
      RETURN
  340 CONTINUE
      IERR = 2
      RETURN
      END
