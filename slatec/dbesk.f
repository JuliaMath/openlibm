*DECK DBESK
      SUBROUTINE DBESK (X, FNU, KODE, N, Y, NZ)
C***BEGIN PROLOGUE  DBESK
C***PURPOSE  Implement forward recursion on the three term recursion
C            relation for a sequence of non-negative order Bessel
C            functions K/SUB(FNU+I-1)/(X), or scaled Bessel functions
C            EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N for real, positive
C            X and non-negative orders FNU.
C***LIBRARY   SLATEC
C***CATEGORY  C10B3
C***TYPE      DOUBLE PRECISION (BESK-S, DBESK-D)
C***KEYWORDS  K BESSEL FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract  **** a double precision routine ****
C         DBESK implements forward recursion on the three term
C         recursion relation for a sequence of non-negative order Bessel
C         functions K/sub(FNU+I-1)/(X), or scaled Bessel functions
C         EXP(X)*K/sub(FNU+I-1)/(X), I=1,..,N for real X .GT. 0.0D0 and
C         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and
C         FNU+1 are obtained from DBSKNU to start the recursion.  If
C         FNU .GE. NULIM, the uniform asymptotic expansion is used for
C         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or
C         70 depending on whether N=1 or N .GE. 2.  Under and overflow
C         tests are made on the leading term of the asymptotic expansion
C         before any extensive computation is done.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         double precision arithmetic.
C
C     Description of Arguments
C
C         Input      X,FNU are double precision
C           X      - X .GT. 0.0D0
C           FNU    - order of the initial K function, FNU .GE. 0.0D0
C           KODE   - a parameter to indicate the scaling option
C                    KODE=1 returns Y(I)=       K/sub(FNU+I-1)/(X),
C                                        I=1,...,N
C                    KODE=2 returns Y(I)=EXP(X)*K/sub(FNU+I-1)/(X),
C                                        I=1,...,N
C           N      - number of members in the sequence, N .GE. 1
C
C         Output     Y is double precision
C           Y      - a vector whose first N components contain values
C                    for the sequence
C                    Y(I)=       k/sub(FNU+I-1)/(X), I=1,...,N  or
C                    Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N
C                    depending on KODE
C           NZ     - number of components of Y set to zero due to
C                    underflow with KODE=1,
C                    NZ=0   , normal return, computation completed
C                    NZ .NE. 0, first NZ components of Y set to zero
C                             due to underflow, Y(I)=0.0D0, I=1,...,NZ
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Overflow - a fatal error
C         Underflow with KODE=1 -  a non-fatal error (NZ .NE. 0)
C
C***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate
C                 or Large Orders, NPL Mathematical Tables 6, Her
C                 Majesty's Stationery Office, London, 1962.
C               N. M. Temme, On the numerical evaluation of the modified
C                 Bessel function of the third kind, Journal of
C                 Computational Physics 19, (1975), pp. 324-337.
C***ROUTINES CALLED  D1MACH, DASYIK, DBESK0, DBESK1, DBSK0E, DBSK1E,
C                    DBSKNU, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790201  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBESK
C
      INTEGER I, J, K, KODE, MZ, N, NB, ND, NN, NUD, NULIM, NZ
      INTEGER I1MACH
      DOUBLE PRECISION CN,DNU,ELIM,ETX,FLGIK,FN,FNN,FNU,GLN,GNU,RTZ,
     1 S, S1, S2, T, TM, TRX, W, X, XLIM, Y, ZN
      DOUBLE PRECISION DBESK0, DBESK1, DBSK1E, DBSK0E, D1MACH
      DIMENSION W(2), NULIM(2), Y(*)
      SAVE NULIM
      DATA NULIM(1),NULIM(2) / 35 , 70 /
C***FIRST EXECUTABLE STATEMENT  DBESK
      NN = -I1MACH(15)
      ELIM = 2.303D0*(NN*D1MACH(5)-3.0D0)
      XLIM = D1MACH(1)*1.0D+3
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 280
      IF (FNU.LT.0.0D0) GO TO 290
      IF (X.LE.0.0D0) GO TO 300
      IF (X.LT.XLIM) GO TO 320
      IF (N.LT.1) GO TO 310
      ETX = KODE - 1
C
C     ND IS A DUMMY VARIABLE FOR N
C     GNU IS A DUMMY VARIABLE FOR FNU
C     NZ = NUMBER OF UNDERFLOWS ON KODE=1
C
      ND = N
      NZ = 0
      NUD = INT(FNU)
      DNU = FNU - NUD
      GNU = FNU
      NN = MIN(2,ND)
      FN = FNU + N - 1
      FNN = FN
      IF (FN.LT.2.0D0) GO TO 150
C
C     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
C     FOR THE LAST ORDER, FNU+N-1.GE.NULIM
C
      ZN = X/FN
      IF (ZN.EQ.0.0D0) GO TO 320
      RTZ = SQRT(1.0D0+ZN*ZN)
      GLN = LOG((1.0D0+RTZ)/ZN)
      T = RTZ*(1.0D0-ETX) + ETX/(ZN+RTZ)
      CN = -FN*(T-GLN)
      IF (CN.GT.ELIM) GO TO 320
      IF (NUD.LT.NULIM(NN)) GO TO 30
      IF (NN.EQ.1) GO TO 20
   10 CONTINUE
C
C     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
C     FOR THE FIRST ORDER, FNU.GE.NULIM
C
      FN = GNU
      ZN = X/FN
      RTZ = SQRT(1.0D0+ZN*ZN)
      GLN = LOG((1.0D0+RTZ)/ZN)
      T = RTZ*(1.0D0-ETX) + ETX/(ZN+RTZ)
      CN = -FN*(T-GLN)
   20 CONTINUE
      IF (CN.LT.-ELIM) GO TO 230
C
C     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM
C
      FLGIK = -1.0D0
      CALL DASYIK(X,GNU,KODE,FLGIK,RTZ,CN,NN,Y)
      IF (NN.EQ.1) GO TO 240
      TRX = 2.0D0/X
      TM = (GNU+GNU+2.0D0)/X
      GO TO 130
C
   30 CONTINUE
      IF (KODE.EQ.2) GO TO 40
C
C     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X)
C     FOR ORDER DNU
C
      IF (X.GT.ELIM) GO TO 230
   40 CONTINUE
      IF (DNU.NE.0.0D0) GO TO 80
      IF (KODE.EQ.2) GO TO 50
      S1 = DBESK0(X)
      GO TO 60
   50 S1 = DBSK0E(X)
   60 CONTINUE
      IF (NUD.EQ.0 .AND. ND.EQ.1) GO TO 120
      IF (KODE.EQ.2) GO TO 70
      S2 = DBESK1(X)
      GO TO 90
   70 S2 = DBSK1E(X)
      GO TO 90
   80 CONTINUE
      NB = 2
      IF (NUD.EQ.0 .AND. ND.EQ.1) NB = 1
      CALL DBSKNU(X, DNU, KODE, NB, W, NZ)
      S1 = W(1)
      IF (NB.EQ.1) GO TO 120
      S2 = W(2)
   90 CONTINUE
      TRX = 2.0D0/X
      TM = (DNU+DNU+2.0D0)/X
C     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
      IF (ND.EQ.1) NUD = NUD - 1
      IF (NUD.GT.0) GO TO 100
      IF (ND.GT.1) GO TO 120
      S1 = S2
      GO TO 120
  100 CONTINUE
      DO 110 I=1,NUD
        S = S2
        S2 = TM*S2 + S1
        S1 = S
        TM = TM + TRX
  110 CONTINUE
      IF (ND.EQ.1) S1 = S2
  120 CONTINUE
      Y(1) = S1
      IF (ND.EQ.1) GO TO 240
      Y(2) = S2
  130 CONTINUE
      IF (ND.EQ.2) GO TO 240
C     FORWARD RECUR FROM FNU+2 TO FNU+N-1
      DO 140 I=3,ND
        Y(I) = TM*Y(I-1) + Y(I-2)
        TM = TM + TRX
  140 CONTINUE
      GO TO 240
C
  150 CONTINUE
C     UNDERFLOW TEST FOR KODE=1
      IF (KODE.EQ.2) GO TO 160
      IF (X.GT.ELIM) GO TO 230
  160 CONTINUE
C     OVERFLOW TEST
      IF (FN.LE.1.0D0) GO TO 170
      IF (-FN*(LOG(X)-0.693D0).GT.ELIM) GO TO 320
  170 CONTINUE
      IF (DNU.EQ.0.0D0) GO TO 180
      CALL DBSKNU(X, FNU, KODE, ND, Y, MZ)
      GO TO 240
  180 CONTINUE
      J = NUD
      IF (J.EQ.1) GO TO 210
      J = J + 1
      IF (KODE.EQ.2) GO TO 190
      Y(J) = DBESK0(X)
      GO TO 200
  190 Y(J) = DBSK0E(X)
  200 IF (ND.EQ.1) GO TO 240
      J = J + 1
  210 IF (KODE.EQ.2) GO TO 220
      Y(J) = DBESK1(X)
      GO TO 240
  220 Y(J) = DBSK1E(X)
      GO TO 240
C
C     UPDATE PARAMETERS ON UNDERFLOW
C
  230 CONTINUE
      NUD = NUD + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 240
      NN = MIN(2,ND)
      GNU = GNU + 1.0D0
      IF (FNN.LT.2.0D0) GO TO 230
      IF (NUD.LT.NULIM(NN)) GO TO 230
      GO TO 10
  240 CONTINUE
      NZ = N - ND
      IF (NZ.EQ.0) RETURN
      IF (ND.EQ.0) GO TO 260
      DO 250 I=1,ND
        J = N - I + 1
        K = ND - I + 1
        Y(J) = Y(K)
  250 CONTINUE
  260 CONTINUE
      DO 270 I=1,NZ
        Y(I) = 0.0D0
  270 CONTINUE
      RETURN
C
C
C
  280 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK',
     +   'SCALING OPTION, KODE, NOT 1 OR 2', 2, 1)
      RETURN
  290 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'ORDER, FNU, LESS THAN ZERO', 2,
     +   1)
      RETURN
  300 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'X LESS THAN OR EQUAL TO ZERO',
     +   2, 1)
      RETURN
  310 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'N LESS THAN ONE', 2, 1)
      RETURN
  320 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK',
     +   'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL', 6, 1)
      RETURN
      END
