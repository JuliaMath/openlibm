*DECK DBSKNU
      SUBROUTINE DBSKNU (X, FNU, KODE, N, Y, NZ)
C***BEGIN PROLOGUE  DBSKNU
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DBESK
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (BESKNU-S, DBSKNU-D)
C***AUTHOR  Amos, D. E., (SNLA)
C***DESCRIPTION
C
C     Abstract  **** A DOUBLE PRECISION routine ****
C         DBSKNU computes N member sequences of K Bessel functions
C         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
C         positive X. Equations of the references are implemented on
C         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X).
C         Forward recursion with the three term recursion relation
C         generates higher orders FNU+I-1, I=1,...,N. The parameter
C         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values
C         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned.
C
C         To start the recursion FNU is normalized to the interval
C         -0.5.LE.DNU.LT.0.5. A special form of the power series is
C         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the
C         K Bessel function in terms of the confluent hypergeometric
C         function U(FNU+0.5,2*FNU+1,X) is implemented on X1.LT.X.LE.X2.
C         For X.GT.X2, the asymptotic expansion for large X is used.
C         When FNU is a half odd integer, a special formula for
C         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
C
C         The maximum number of significant digits obtainable
C         is the smaller of 14 and the number of digits carried in
C         DOUBLE PRECISION arithmetic.
C
C         DBSKNU assumes that a significant digit SINH function is
C         available.
C
C     Description of Arguments
C
C         INPUT      X,FNU are DOUBLE PRECISION
C           X      - X.GT.0.0D0
C           FNU    - Order of initial K function, FNU.GE.0.0D0
C           N      - Number of members of the sequence, N.GE.1
C           KODE   - A parameter to indicate the scaling option
C                    KODE= 1  returns
C                             Y(I)=       K/SUB(FNU+I-1)/(X)
C                                  I=1,...,N
C                        = 2  returns
C                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X)
C                                  I=1,...,N
C
C         OUTPUT     Y is DOUBLE PRECISION
C           Y      - A vector whose first N components contain values
C                    for the sequence
C                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or
C                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N
C                    depending on KODE
C           NZ     - Number of components set to zero due to
C                    underflow,
C                    NZ= 0   , normal return
C                    NZ.NE.0 , first NZ components of Y set to zero
C                              due to underflow, Y(I)=0.0D0,I=1,...,NZ
C
C     Error Conditions
C         Improper input arguments - a fatal error
C         Overflow - a fatal error
C         Underflow with KODE=1 - a non-fatal error (NZ.NE.0)
C
C***SEE ALSO  DBESK
C***REFERENCES  N. M. Temme, On the numerical evaluation of the modified
C                 Bessel function of the third kind, Journal of
C                 Computational Physics 19, (1975), pp. 324-337.
C***ROUTINES CALLED  D1MACH, DGAMMA, I1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790201  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C   900727  Added EXTERNAL statement.  (WRB)
C   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DBSKNU
C
      INTEGER I, IFLAG, INU, J, K, KK, KODE, KODED, N, NN, NZ
      INTEGER I1MACH
      DOUBLE PRECISION A,AK,A1,A2,B,BK,CC,CK,COEF,CX,DK,DNU,DNU2,ELIM,
     1 ETEST, EX, F, FC, FHS, FK, FKS, FLRX, FMU, FNU, G1, G2, P, PI,
     2 PT, P1, P2, Q, RTHPI, RX, S, SMU, SQK, ST, S1, S2, TM, TOL, T1,
     3 T2, X, X1, X2, Y
      DIMENSION A(160), B(160), Y(*), CC(8)
      DOUBLE PRECISION DGAMMA, D1MACH
      EXTERNAL DGAMMA
      SAVE X1, X2, PI, RTHPI, CC
      DATA X1, X2 / 2.0D0, 17.0D0 /
      DATA PI,RTHPI        / 3.14159265358979D+00, 1.25331413731550D+00/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)
     1                     / 5.77215664901533D-01,-4.20026350340952D-02,
     2-4.21977345555443D-02, 7.21894324666300D-03,-2.15241674114900D-04,
     3-2.01348547807000D-05, 1.13302723200000D-06, 6.11609500000000D-09/
C***FIRST EXECUTABLE STATEMENT  DBSKNU
      KK = -I1MACH(15)
      ELIM = 2.303D0*(KK*D1MACH(5)-3.0D0)
      AK = D1MACH(3)
      TOL = MAX(AK,1.0D-15)
      IF (X.LE.0.0D0) GO TO 350
      IF (FNU.LT.0.0D0) GO TO 360
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 370
      IF (N.LT.1) GO TO 380
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RX = 2.0D0/X
      INU = INT(FNU+0.5D0)
      DNU = FNU - INU
      IF (ABS(DNU).EQ.0.5D0) GO TO 120
      DNU2 = 0.0D0
      IF (ABS(DNU).LT.TOL) GO TO 10
      DNU2 = DNU*DNU
   10 CONTINUE
      IF (X.GT.X1) GO TO 120
C
C     SERIES FOR X.LE.X1
C
      A1 = 1.0D0 - DNU
      A2 = 1.0D0 + DNU
      T1 = 1.0D0/DGAMMA(A1)
      T2 = 1.0D0/DGAMMA(A2)
      IF (ABS(DNU).GT.0.1D0) GO TO 40
C     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
      S = CC(1)
      AK = 1.0D0
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (ABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5D0
      SMU = 1.0D0
      FC = 1.0D0
      FLRX = LOG(RX)
      FMU = DNU*FLRX
      IF (DNU.EQ.0.0D0) GO TO 60
      FC = DNU*PI
      FC = FC/SIN(FC)
      IF (FMU.NE.0.0D0) SMU = SINH(FMU)/FMU
   60 CONTINUE
      F = FC*(G1*COSH(FMU)+G2*FLRX*SMU)
      FC = EXP(FMU)
      P = 0.5D0*FC/T2
      Q = 0.5D0/(FC*T1)
      AK = 1.0D0
      CK = 1.0D0
      BK = 1.0D0
      S1 = F
      S2 = P
      IF (INU.GT.0 .OR. N.GT.1) GO TO 90
      IF (X.LT.TOL) GO TO 80
      CX = X*X*0.25D0
   70 CONTINUE
      F = (AK*F+P+Q)/(BK-DNU2)
      P = P/(AK-DNU)
      Q = Q/(AK+DNU)
      CK = CK*CX/AK
      T1 = CK*F
      S1 = S1 + T1
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      S = ABS(T1)/(1.0D0+ABS(S1))
      IF (S.GT.TOL) GO TO 70
   80 CONTINUE
      Y(1) = S1
      IF (KODED.EQ.1) RETURN
      Y(1) = S1*EXP(X)
      RETURN
   90 CONTINUE
      IF (X.LT.TOL) GO TO 110
      CX = X*X*0.25D0
  100 CONTINUE
      F = (AK*F+P+Q)/(BK-DNU2)
      P = P/(AK-DNU)
      Q = Q/(AK+DNU)
      CK = CK*CX/AK
      T1 = CK*F
      S1 = S1 + T1
      T2 = CK*(P-AK*F)
      S2 = S2 + T2
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      S = ABS(T1)/(1.0D0+ABS(S1)) + ABS(T2)/(1.0D0+ABS(S2))
      IF (S.GT.TOL) GO TO 100
  110 CONTINUE
      S2 = S2*RX
      IF (KODED.EQ.1) GO TO 170
      F = EXP(X)
      S1 = S1*F
      S2 = S2*F
      GO TO 170
  120 CONTINUE
      COEF = RTHPI/SQRT(X)
      IF (KODED.EQ.2) GO TO 130
      IF (X.GT.ELIM) GO TO 330
      COEF = COEF*EXP(-X)
  130 CONTINUE
      IF (ABS(DNU).EQ.0.5D0) GO TO 340
      IF (X.GT.X2) GO TO 280
C
C     MILLER ALGORITHM FOR X1.LT.X.LE.X2
C
      ETEST = COS(PI*DNU)/(PI*X*TOL)
      FKS = 1.0D0
      FHS = 0.25D0
      FK = 0.0D0
      CK = X + X + 2.0D0
      P1 = 0.0D0
      P2 = 1.0D0
      K = 0
  140 CONTINUE
      K = K + 1
      FK = FK + 1.0D0
      AK = (FHS-DNU2)/(FKS+FK)
      BK = CK/(FK+1.0D0)
      PT = P2
      P2 = BK*P2 - AK*P1
      P1 = PT
      A(K) = AK
      B(K) = BK
      CK = CK + 2.0D0
      FKS = FKS + FK + FK + 1.0D0
      FHS = FHS + FK + FK
      IF (ETEST.GT.FK*P1) GO TO 140
      KK = K
      S = 1.0D0
      P1 = 0.0D0
      P2 = 1.0D0
      DO 150 I=1,K
        PT = P2
        P2 = (B(KK)*P2-P1)/A(KK)
        P1 = PT
        S = S + P2
        KK = KK - 1
  150 CONTINUE
      S1 = COEF*(P2/S)
      IF (INU.GT.0 .OR. N.GT.1) GO TO 160
      GO TO 200
  160 CONTINUE
      S2 = S1*(X+DNU+0.5D0-P1/P2)/X
C
C     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
C
  170 CONTINUE
      CK = (DNU+DNU+2.0D0)/X
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 180
      IF (N.GT.1) GO TO 200
      S1 = S2
      GO TO 200
  180 CONTINUE
      DO 190 I=1,INU
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        CK = CK + RX
  190 CONTINUE
      IF (N.EQ.1) S1 = S2
  200 CONTINUE
      IF (IFLAG.EQ.1) GO TO 220
      Y(1) = S1
      IF (N.EQ.1) RETURN
      Y(2) = S2
      IF (N.EQ.2) RETURN
      DO 210 I=3,N
        Y(I) = CK*Y(I-1) + Y(I-2)
        CK = CK + RX
  210 CONTINUE
      RETURN
C     IFLAG=1 CASES
  220 CONTINUE
      S = -X + LOG(S1)
      Y(1) = 0.0D0
      NZ = 1
      IF (S.LT.-ELIM) GO TO 230
      Y(1) = EXP(S)
      NZ = 0
  230 CONTINUE
      IF (N.EQ.1) RETURN
      S = -X + LOG(S2)
      Y(2) = 0.0D0
      NZ = NZ + 1
      IF (S.LT.-ELIM) GO TO 240
      NZ = NZ - 1
      Y(2) = EXP(S)
  240 CONTINUE
      IF (N.EQ.2) RETURN
      KK = 2
      IF (NZ.LT.2) GO TO 260
      DO 250 I=3,N
        KK = I
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        CK = CK + RX
        S = -X + LOG(S2)
        NZ = NZ + 1
        Y(I) = 0.0D0
        IF (S.LT.-ELIM) GO TO 250
        Y(I) = EXP(S)
        NZ = NZ - 1
        GO TO 260
  250 CONTINUE
      RETURN
  260 CONTINUE
      IF (KK.EQ.N) RETURN
      S2 = S2*CK + S1
      CK = CK + RX
      KK = KK + 1
      Y(KK) = EXP(-X+LOG(S2))
      IF (KK.EQ.N) RETURN
      KK = KK + 1
      DO 270 I=KK,N
        Y(I) = CK*Y(I-1) + Y(I-2)
        CK = CK + RX
  270 CONTINUE
      RETURN
C
C     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2
C
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
  280 CONTINUE
      NN = 2
      IF (INU.EQ.0 .AND. N.EQ.1) NN = 1
      DNU2 = DNU + DNU
      FMU = 0.0D0
      IF (ABS(DNU2).LT.TOL) GO TO 290
      FMU = DNU2*DNU2
  290 CONTINUE
      EX = X*8.0D0
      S2 = 0.0D0
      DO 320 K=1,NN
        S1 = S2
        S = 1.0D0
        AK = 0.0D0
        CK = 1.0D0
        SQK = 1.0D0
        DK = EX
        DO 300 J=1,30
          CK = CK*(FMU-SQK)/DK
          S = S + CK
          DK = DK + EX
          AK = AK + 8.0D0
          SQK = SQK + AK
          IF (ABS(CK).LT.TOL) GO TO 310
  300   CONTINUE
  310   S2 = S*COEF
        FMU = FMU + 8.0D0*DNU + 4.0D0
  320 CONTINUE
      IF (NN.GT.1) GO TO 170
      S1 = S2
      GO TO 200
  330 CONTINUE
      KODED = 2
      IFLAG = 1
      GO TO 120
C
C     FNU=HALF ODD INTEGER CASE
C
  340 CONTINUE
      S1 = COEF
      S2 = COEF
      GO TO 170
C
C
  350 CALL XERMSG ('SLATEC', 'DBSKNU', 'X NOT GREATER THAN ZERO', 2, 1)
      RETURN
  360 CALL XERMSG ('SLATEC', 'DBSKNU', 'FNU NOT ZERO OR POSITIVE', 2,
     +   1)
      RETURN
  370 CALL XERMSG ('SLATEC', 'DBSKNU', 'KODE NOT 1 OR 2', 2, 1)
      RETURN
  380 CALL XERMSG ('SLATEC', 'DBSKNU', 'N NOT GREATER THAN 0', 2, 1)
      RETURN
      END
