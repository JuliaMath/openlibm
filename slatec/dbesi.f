*DECK DBESI
      SUBROUTINE DBESI (X, ALPHA, KODE, N, Y, NZ)
C***BEGIN PROLOGUE  DBESI
C***PURPOSE  Compute an N member sequence of I Bessel functions
C            I/SUB(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
C            EXP(-X)*I/SUB(ALPHA+K-1)/(X), K=1,...,N for nonnegative
C            ALPHA and X.
C***LIBRARY   SLATEC
C***CATEGORY  C10B3
C***TYPE      DOUBLE PRECISION (BESI-S, DBESI-D)
C***KEYWORDS  I BESSEL FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNLA)
C           Daniel, S. L., (SNLA)
C***DESCRIPTION
C
C     Abstract  **** a double precision routine ****
C         DBESI computes an N member sequence of I Bessel functions
C         I/sub(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
C         EXP(-X)*I/sub(ALPHA+K-1)/(X), K=1,...,N for nonnegative ALPHA
C         and X.  A combination of the power series, the asymptotic
C         expansion for X to infinity, and the uniform asymptotic
C         expansion for NU to infinity are applied over subdivisions of
C         the (NU,X) plane.  For values not covered by one of these
C         formulae, the order is incremented by an integer so that one
C         of these formulae apply.  Backward recursion is used to reduce
C         orders by integer values.  The asymptotic expansion for X to
C         infinity is used only when the entire sequence (specifically
C         the last member) lies within the region covered by the
C         expansion.  Leading terms of these expansions are used to test
C         for over or underflow where appropriate.  If a sequence is
C         requested and the last member would underflow, the result is
C         set to zero and the next lower order tried, etc., until a
C         member comes on scale or all are set to zero.  An overflow
C         cannot occur with scaling.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C     Description of Arguments
C
C         Input      X,ALPHA are double precision
C           X      - X .GE. 0.0D0
C           ALPHA  - order of first member of the sequence,
C                    ALPHA .GE. 0.0D0
C           KODE   - a parameter to indicate the scaling option
C                    KODE=1 returns
C                           Y(K)=        I/sub(ALPHA+K-1)/(X),
C                                K=1,...,N
C                    KODE=2 returns
C                           Y(K)=EXP(-X)*I/sub(ALPHA+K-1)/(X),
C                                K=1,...,N
C           N      - number of members in the sequence, N .GE. 1
C
C         Output     Y is double precision
C           Y      - a vector whose first N components contain
C                    values for I/sub(ALPHA+K-1)/(X) or scaled
C                    values for EXP(-X)*I/sub(ALPHA+K-1)/(X),
C                    K=1,...,N depending on KODE
C           NZ     - number of components of Y set to zero due to
C                    underflow,
C                    NZ=0   , normal return, computation completed
C                    NZ .NE. 0, last NZ components of Y set to zero,
C                             Y(K)=0.0D0, K=N-NZ+1,...,N.
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Overflow with KODE=1 - a fatal error
C         Underflow - a non-fatal error(NZ .NE. 0)
C
C***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
C                 subroutines IBESS and JBESS for Bessel functions
C                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
C                 Transactions on Mathematical Software 3, (1977),
C                 pp. 76-92.
C               F. W. J. Olver, Tables of Bessel Functions of Moderate
C                 or Large Orders, NPL Mathematical Tables 6, Her
C                 Majesty's Stationery Office, London, 1962.
C***ROUTINES CALLED  D1MACH, DASYIK, DLNGAM, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBESI
C
      INTEGER I, IALP, IN, INLIM, IS, I1, K, KK, KM, KODE, KT,
     1 N, NN, NS, NZ
      INTEGER I1MACH
      DOUBLE PRECISION AIN,AK,AKM,ALPHA,ANS,AP,ARG,ATOL,TOLLN,DFN,
     1 DTM, DX, EARG, ELIM, ETX, FLGIK,FN, FNF, FNI,FNP1,FNU,GLN,RA,
     2 RTTPI, S, SX, SXO2, S1, S2, T, TA, TB, TEMP, TFN, TM, TOL,
     3 TRX, T2, X, XO2, XO2L, Y, Z
      DOUBLE PRECISION D1MACH, DLNGAM
      DIMENSION Y(*), TEMP(3)
      SAVE RTTPI, INLIM
      DATA RTTPI           / 3.98942280401433D-01/
      DATA INLIM           /          80         /
C***FIRST EXECUTABLE STATEMENT  DBESI
      NZ = 0
      KT = 1
C     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
C     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
      RA = D1MACH(3)
      TOL = MAX(RA,1.0D-15)
      I1 = -I1MACH(15)
      GLN = D1MACH(5)
      ELIM = 2.303D0*(I1*GLN-3.0D0)
C     TOLLN = -LN(TOL)
      I1 = I1MACH(14)+1
      TOLLN = 2.303D0*GLN*I1
      TOLLN = MIN(TOLLN,34.5388D0)
      IF (N-1) 590, 10, 20
   10 KT = 2
   20 NN = N
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 570
      IF (X) 600, 30, 80
   30 IF (ALPHA) 580, 40, 50
   40 Y(1) = 1.0D0
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.0D0
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.0D0) GO TO 580
C
      IALP = INT(ALPHA)
      FNI = IALP + N - 1
      FNF = ALPHA - IALP
      DFN = FNI + FNF
      FNU = DFN
      IN = 0
      XO2 = X*0.5D0
      SXO2 = XO2*XO2
      ETX = KODE - 1
      SX = ETX*X
C
C     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
C     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
C     APPLIED.
C
      IF (SXO2.LE.(FNU+1.0D0)) GO TO 90
      IF (X.LE.12.0D0) GO TO 110
      FN = 0.55D0*FNU*FNU
      FN = MAX(17.0D0,FN)
      IF (X.GE.FN) GO TO 430
      ANS = MAX(36.0D0-FNU,0.0D0)
      NS = INT(ANS)
      FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      KM = N - 1 + NS
      IF (KM.GT.0) IS = 3
      GO TO 120
   90 FN = FNU
      FNP1 = FN + 1.0D0
      XO2L = LOG(XO2)
      IS = KT
      IF (X.LE.0.5D0) GO TO 230
      NS = 0
  100 FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.0D0
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 230
  110 XO2L = LOG(XO2)
      NS = INT(SXO2-FNU)
      GO TO 100
  120 CONTINUE
C
C     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
C
      IF (KODE.EQ.2) GO TO 130
      IF (ALPHA.LT.1.0D0) GO TO 150
      Z = X/ALPHA
      RA = SQRT(1.0D0+Z*Z)
      GLN = LOG((1.0D0+RA)/Z)
      T = RA*(1.0D0-ETX) + ETX/(Z+RA)
      ARG = ALPHA*(T-GLN)
      IF (ARG.GT.ELIM) GO TO 610
      IF (KM.EQ.0) GO TO 140
  130 CONTINUE
C
C     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
C
      Z = X/FN
      RA = SQRT(1.0D0+Z*Z)
      GLN = LOG((1.0D0+RA)/Z)
      T = RA*(1.0D0-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  140 IF (ARG.LT.(-ELIM)) GO TO 280
      GO TO 190
  150 IF (X.GT.ELIM) GO TO 610
      GO TO 130
C
C     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
C
  160 IF (KM.NE.0) GO TO 170
      Y(1) = TEMP(3)
      RETURN
  170 TEMP(1) = TEMP(3)
      IN = NS
      KT = 1
      I1 = 0
  180 CONTINUE
      IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 350
      Z = X/FN
      RA = SQRT(1.0D0+Z*Z)
      GLN = LOG((1.0D0+RA)/Z)
      T = RA*(1.0D0-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  190 CONTINUE
      I1 = ABS(3-IS)
      I1 = MAX(I1,1)
      FLGIK = 1.0D0
      CALL DASYIK(X,FN,KODE,FLGIK,RA,ARG,I1,TEMP(IS))
      GO TO (180, 350, 510), IS
C
C     SERIES FOR (X/2)**2.LE.NU+1
C
  230 CONTINUE
      GLN = DLNGAM(FNP1)
      ARG = FN*XO2L - GLN - SX
      IF (ARG.LT.(-ELIM)) GO TO 300
      EARG = EXP(ARG)
  240 CONTINUE
      S = 1.0D0
      IF (X.LT.TOL) GO TO 260
      AK = 3.0D0
      T2 = 1.0D0
      T = 1.0D0
      S1 = FN
      DO 250 K=1,17
        S2 = T2 + S1
        T = T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 260
        T2 = T2 + AK
        AK = AK + 2.0D0
        S1 = S1 + FN
  250 CONTINUE
  260 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (270, 350, 500), IS
  270 EARG = EARG*FN/XO2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 240
C
C     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
C
  280 Y(NN) = 0.0D0
      NN = NN - 1
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 340, 290, 130
  290 KT = 2
      IS = 2
      GO TO 130
  300 Y(NN) = 0.0D0
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF (NN-1) 340, 310, 320
  310 KT = 2
      IS = 2
  320 IF (SXO2.LE.FNP1) GO TO 330
      GO TO 130
  330 ARG = ARG - XO2L + LOG(FNP1)
      IF (ARG.LT.(-ELIM)) GO TO 300
      GO TO 230
  340 NZ = N - NN
      RETURN
C
C     BACKWARD RECURSION SECTION
C
  350 CONTINUE
      NZ = N - NN
  360 CONTINUE
      IF(KT.EQ.2) GO TO 420
      S1 = TEMP(1)
      S2 = TEMP(2)
      TRX = 2.0D0/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      IF (IN.EQ.0) GO TO 390
C     BACKWARD RECUR TO INDEX ALPHA+NN-1
      DO 380 I=1,IN
        S = S2
        S2 = TM*S2 + S1
        S1 = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  380 CONTINUE
      Y(NN) = S1
      IF (NN.EQ.1) RETURN
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
      GO TO 400
  390 CONTINUE
C     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = S1
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
  400 K = NN + 1
      DO 410 I=3,NN
        K = K - 1
        Y(K-2) = TM*Y(K-1) + Y(K)
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  410 CONTINUE
      RETURN
  420 Y(1) = TEMP(2)
      RETURN
C
C     ASYMPTOTIC EXPANSION FOR X TO INFINITY
C
  430 CONTINUE
      EARG = RTTPI/SQRT(X)
      IF (KODE.EQ.2) GO TO 440
      IF (X.GT.ELIM) GO TO 610
      EARG = EARG*EXP(X)
  440 ETX = 8.0D0*X
      IS = KT
      IN = 0
      FN = FNU
  450 DX = FNI + FNI
      TM = 0.0D0
      IF (FNI.EQ.0.0D0 .AND. ABS(FNF).LT.TOL) GO TO 460
      TM = 4.0D0*FNF*(FNI+FNI+FNF)
  460 CONTINUE
      DTM = DX*DX
      S1 = ETX
      TRX = DTM - 1.0D0
      DX = -(TRX+TM)/ETX
      T = DX
      S = 1.0D0 + DX
      ATOL = TOL*ABS(S)
      S2 = 1.0D0
      AK = 8.0D0
      DO 470 K=1,25
        S1 = S1 + ETX
        S2 = S2 + AK
        DX = DTM - S2
        AP = DX + TM
        T = -T*AP/S1
        S = S + T
        IF (ABS(T).LE.ATOL) GO TO 480
        AK = AK + 8.0D0
  470 CONTINUE
  480 TEMP(IS) = S*EARG
      IF(IS.EQ.2) GO TO 360
      IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      GO TO 450
C
C     BACKWARD RECURSION WITH NORMALIZATION BY
C     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
C
  500 CONTINUE
C     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = MAX(3.0D0-FN,0.0D0)
      KM = INT(AKM)
      TFN = FN + KM
      TA = (GLN+TFN-0.9189385332D0-0.0833333333D0/TFN)/(TFN+0.5D0)
      TA = XO2L - TA
      TB = -(1.0D0-1.0D0/TFN)/TFN
      AIN = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.5D0
      IN = INT(AIN)
      IN = IN + KM
      GO TO 520
  510 CONTINUE
C     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      T = 1.0D0/(FN*RA)
      AIN = TOLLN/(GLN+SQRT(GLN*GLN+T*TOLLN)) + 1.5D0
      IN = INT(AIN)
      IF (IN.GT.INLIM) GO TO 160
  520 CONTINUE
      TRX = 2.0D0/X
      DTM = FNI + IN
      TM = (DTM+FNF)*TRX
      TA = 0.0D0
      TB = TOL
      KK = 1
  530 CONTINUE
C
C     BACKWARD RECUR UNINDEXED
C
      DO 540 I=1,IN
        S = TB
        TB = TM*TB + TA
        TA = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  540 CONTINUE
C     NORMALIZATION
      IF (KK.NE.1) GO TO 550
      TA = (TA/TB)*TEMP(3)
      TB = TEMP(3)
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 530
  550 Y(NN) = TB
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      TB = TM*TB + TA
      K = NN - 1
      Y(K) = TB
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.0D0
      TM = (DTM+FNF)*TRX
      KM = K - 1
C
C     BACKWARD RECUR INDEXED
C
      DO 560 I=1,KM
        Y(K-1) = TM*Y(K) + Y(K+1)
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  560 CONTINUE
      RETURN
C
C
C
  570 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI',
     +   'SCALING OPTION, KODE, NOT 1 OR 2.', 2, 1)
      RETURN
  580 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'ORDER, ALPHA, LESS THAN ZERO.',
     +   2, 1)
      RETURN
  590 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'N LESS THAN ONE.', 2, 1)
      RETURN
  600 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'X LESS THAN ZERO.', 2, 1)
      RETURN
  610 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI',
     +   'OVERFLOW, X TOO LARGE FOR KODE = 1.', 6, 1)
      RETURN
      END
